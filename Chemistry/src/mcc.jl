using ParticleInCell
module MCC
	using ParticleInCell
	struct ElasticIsotropic end
	struct ElasticBackward end
	struct Ionization energy :: Float64 end
	struct Excitation energy :: Float64 end
	struct Collision{T}
		type :: T
		rate
		source :: KineticSpecies
		target :: FluidSpecies
		products :: Vector{KineticSpecies}
	end
end

mutable struct MonteCarloCollisions
	collisions :: Vector{MCC.Collision}
	max_σₜg :: Float64
	m :: Float64 # mass in eV
	remainder :: Float64
	#ε :: Vector{Float64}
end

@inline mass(species) = species.m/1.60217646e-19 # qe = 1.60217646e-19C
function MonteCarloCollisions(collisions :: Vector{MCC.Collision})
	ε = [0.0]
	for collision in collisions
		εᵢ = collision.rate.nodes[:,1]
		println(collision, ": ", extrema(εᵢ))
		append!(ε, εᵢ)
	end
	ε = sort(unique(ε))

	collision = first(collisions)

	mₛ = mass(collision.source)
	α = sqrt(2.0/mₛ)
	v = α .* sqrt.(ε)

	σₜg = zero(v)
	for collision in collisions
		σᵢ = collision.rate
		σₜg .+= σᵢ.(ε) .* v
	end
	println(collision.source, " MCC Energy range: ", extrema(ε))
	println(collision.source, " MCC Velocity range: ", extrema(v))
	println(collision.source, " MCC maximal velocity: ", v[argmax(σₜg)])
	return MonteCarloCollisions(collisions, maximum(σₜg), mₛ, 0.0)
end

@inline function isotropic_distribution()
    χ = 2π*rand()
	sinχ = sin(χ)
    cosχ = cos(χ)
    
    η = 2π*rand()
    sinη = sin(η)
    cosη = cos(η)
    return sinχ, cosχ, sinη, cosη
end

@inline function cosine_distribution()
	sinχ = sqrt(rand())
    cosχ =-sqrt(1.0 - sinχ^2)
    
    η = 2π*rand()
    sinη = sin(η)
    cosη = cos(η)
    return sinχ, cosχ, sinη, cosη
end

function thermal_speed(T, m)
	kB = 1.3806503e-23
	sqrt(2kB*T/m)
end

function maxwellian_velocity(ν)
    v = [randn(), randn(), 0]
	v ./ norm(v) .* ν
end

function unrotated(sinθ, cosθ, sinϕ, cosϕ)
	[ cosϕ*cosθ -sinϕ*cosθ -sinθ
	 -sinϕ       cosϕ        0.0
	  cosϕ*sinθ  sinϕ*sinθ  cosθ]
end

function euler_angles(v)
	vx = v[1]
	vy = v[2]
	vz = v[3]
	nv = norm(v)

	cosθ = vz/nv
	sinθ = sqrt(1.0 - cosθ^2)
	
	if sinθ ≈ 0.0
		cosϕ = 1.0
		sinϕ = 0.0
	else
		cosϕ = vx/nv/sinθ
		sinϕ = vy/nv/sinθ
	end
	return sinθ, cosθ, sinϕ, cosϕ
end

function isotropic_scattering(v)
	sinθ, cosθ, sinϕ, cosϕ = euler_angles(v)
    sinχ, cosχ, sinη, cosη = isotropic_distribution()
	T = unrotated(sinθ, cosθ, sinϕ, cosϕ)
    vec([sinχ*cosη sinχ*sinη cosχ] * T)
end

function diffuse_reflection(v)
	sinθ, cosθ, sinϕ, cosϕ = euler_angles(v)
	sinχ, cosχ, sinη, cosη = cosine_distribution()
    T = unrotated(sinθ, cosθ, sinϕ, cosϕ)
	vec([sinχ*cosη sinχ*sinη cosχ] * T)
end

function perform!(collision::MCC.Collision{MCC.ElasticBackward}, p, Δt, grid)
	source, target = collision.source, collision.target
	mr1 = source.m/(source.m + target.m)
	mr2 = target.m/(source.m + target.m)
	tv = maxwellian_velocity(thermal_speed(target.T, target.m))
	sv = view(source.v, p, :)

	vr = sv .- tv
	w  = mr1 .* sv .+ mr2 .* tv
	vr_cp = norm(vr) * diffuse_reflection(vr)

	sv .= w .+ mr2 * vr_cp;
	tv .= w .- mr1 * vr_cp;
end

function perform!(collision::MCC.Collision{MCC.ElasticIsotropic}, p, Δt, grid)
	source, target = collision.source, collision.target
	mr1 = source.m/(source.m + target.m)
	mr2 = target.m/(source.m + target.m)
	tv = maxwellian_velocity(thermal_speed(target.T, target.m))
	sv = view(source.v, p, :)

	vr = sv .- tv
	w  = mr1 .* sv .+ mr2 .* tv
	
	vr_cp = norm(vr) * isotropic_scattering(vr)
	sv .= w .+ mr2 * vr_cp;
	tv .= w .- mr1 * vr_cp;
end

function perform!(collision::MCC.Collision{MCC.Ionization}, p, Δt, grid)
	source, target = collision.source, collision.target
	mₛ = mass(source)
	sv = view(source.v, p, :)
	sE = 0.5mₛ*dot(sv, sv) - collision.type.energy
	if sE < 0
	  println("Source species has not enough energy for ionization: ", sE)
	  return
	end

	# randomly redistribute the remaining energy to the two electrons
	e1E = sE * rand()
	e2E = sE - e1E
	# incident electron's speed reduced by the ionization energy
	α = sqrt(2.0/mₛ)
	e1ν = α*sqrt(e1E)
	e2ν = α*sqrt(e2E)

	sv .= e1ν .* diffuse_reflection(sv)

	# create new electron backscattered from the gas atom
	source.x[source.np+1,:] .= source.x[p,:]
	source.v[source.np+1,:] .= e2ν .* diffuse_reflection(sv)
	source.np += 1

	# create new ion with 
	tv = maxwellian_velocity(thermal_speed(target.T, target.m))
	for product in collision.products
		if product == source
			continue
		end
		mp = source.w0/product.w0 + rand()
		# assume the new electron and ion are created at the neutral temperature
		for i=1:round(Integer, mp)
			product.x[product.np+1,:] .= source.x[p,:]
			product.v[product.np+1,:] .= tv
			product.np += 1
		end
	end
end

function perform!(collision::MCC.Collision{MCC.Excitation}, p, Δt, grid)
	source, target = collision.source, collision.target
	mₛ = mass(source)
	sv = view(source.v, p, :)
	sE = 0.5mₛ*dot(sv, sv) - collision.type.energy;
	if sE < 0.0
	  println("Source species has not enough energy for excitation: ", sE)
	  return
	end
	
	# speed reduced by the excitation energy
	α = sqrt(2/mₛ)
	eν = α*sqrt(sE) 
	sv.= eν .* isotropic_scattering(sv)
end

function PIC.perform!(mcc::MonteCarloCollisions, E, Δt, config)
	nx, ny = size(config.grid)
	N = length(mcc.collisions)
	ν = zeros(nx, ny, N) # collision frequency
	Δx, Δy = config.grid.Δh

	collision = first(mcc.collisions)
	density = PIC.density(collision.target, config.grid)
	source, target = collision.source, collision.target
	max_σₜg = mcc.max_σₜg

	max_n₀ = maximum(density)
	max_Pₜ = 1.0 - exp(-max_n₀ * max_σₜg * Δt)
	if max_Pₜ > 1.0/N
		@assert false "Maximum probability ($max_Pₜ) is greater than 1/$N"
	end
	
	mcc.remainder, Nc = modf(N * max_Pₜ * source.np + mcc.remainder)
	#println("Checking ", Nc, " particles (", source, ")")
	#println("...plus ", mcc.remainder, " in the next iteration")
	for ~=1:Nc
		# choose a particle
		p = rand(1:source.np)
		i, j, _, _ = PIC.particle_cell(source.x, p, config.grid.Δh)
		n = density[i,j]
		if n < 0
			println("Density is negative, skipping")
			continue
		end

		# select a collision process
		U = rand()
		k = floor(Int64, N*U + 1)
		collision = mcc.collisions[k]

		# calculate the energy of the particle
		sv = source.v
		tv = (target.q/target.m)*E*Δt
		g  = norm(tv[i,j,:] .- sv[p,:])
		ε  = 0.5mcc.m * g^2
		σₖg = collision.rate(ε) * g

		Pₖ  = 1.0 .- exp.(-n .* σₖg .* Δt)
		Pₖ /= N * max_Pₜ
		if Pₖ > 1.0
			println("Particle's velocity: ", sv[p,:])
			println("Gas thermal velocity: ", tv[i,j,:])
			println("Gas density: ", n)
			println("Energy: ", ε)
			@assert false "Energy outside of the range"
		end

		if U > k/N - Pₖ
			perform!(collision, p, Δt, config.grid)
			ν[i,j,k] += 1
		end
	end
	println("Collisions: ", maximum(sum(ν;dims=3)), " out of ", Nc)
	for k=1:N
		@field "nuMCC-"*source.name*"-"*string(k) "1/m^2" ν[:,:,k] config.grid
	end
end

function accept(reaction::ChemicalReaction)
	default_type = MCC.ElasticIsotropic()
	products = KineticSpecies[]
	source = target = nothing
	@assert length(reaction.reactants) == 2 "Monte Carlo Collisions support only two reacting species: one fluid and one kinetic"
	for (r, _) in reaction.reactants
		if PIC.is_fluid(r) target = r else source = r end
	end
	for (p, c) in reaction.stoichiometry
		if c > 0 push!(products, p) end
	end

	if source == nothing error("Reaction without particle species") end
	if target == nothing error("Reaction without fluid species") end

	if isnothing(reaction.type)
		MCC.Collision(default_type, reaction.rate, source, target, products)
	else
		MCC.Collision(reaction.type, reaction.rate, source, target, products)
	end
end

function mcc(reactions)
	collisions = MCC.Collision[]
	for reaction in reactions
		collision = accept(reaction)
		push!(collisions, collision)
	end
	MonteCarloCollisions(collisions)
end

Base.show(io :: IO, r :: MCC.Collision{MCC.ElasticBackward}) = print(io, r.rate, ", ", r.source, " + ", r.target, "-->", r.source, " + ", r.target, "\telastic (backward)") 
Base.show(io :: IO, r :: MCC.Collision{MCC.ElasticIsotropic}) = print(io, r.rate, ", ", r.source, " + ", r.target, "-->", r.source, " + ", r.target, "\telastic (isotropic)") 
Base.show(io :: IO, r :: MCC.Collision{MCC.Excitation}) = print(io, r.rate, ", ", r.source, " + ", r.target, "-->", r.source, " + ", r.target, "\texcitation") 
Base.show(io :: IO, r :: MCC.Collision{MCC.Ionization}) = print(io, r.rate, ", ", r.source, " + ", r.target, "-->", r.products, "\tionization")