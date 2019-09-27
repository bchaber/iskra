module MCC
	struct ElasticCollision
		rate
		source # particles
		target # fluid
	end

	struct IonizationCollision
		rate
		source # particles
		target # fluid
		products # particles
	end
end

struct MonteCarloCollisions
	collisions
end

function isotropic_velocity(ν)
	θ = 2rand() * π
	r = 2rand() - 1
	a = sqrt(1 - r^2)
	[cos(θ)*a, sin(θ)*a, r] .* ν
end

function thermal_velocity(T, m)
	kB = 1.3806503e-23
	sqrt(2kB*T/m)
end

function maxwellian_velocity(ν)
    v = [randn(), randn(), 0]
	v ./ norm(v) .* ν
end

function perform!(collision::MCC.ElasticCollision, p, Δt, grid)
	ν = norm(collision.source.v[p,:])
	collision.source.v[p,:] .= maxwellian_velocity(ν)
end

function perform!(collision::MCC.IonizationCollision, p, Δt, grid)
	source = collision.source
	qe = 1.60217646e-19
	sv = view(source.v, p, :)
	sE = 0.5source.m*dot(sv, sv)/qe - 12.0697 #target.ionization_energy;
	if sE < 0
	  println("Source species has not enough energy for ionization: ", sE)
	  return
	end
	# randomly redistribute the remaining energy to the two electrons
	e1E = sE * rand()
	e2E = sE - e1E
	# speed reduced by the ionization energy
	e1ν = sqrt(e1E*qe*2/source.m)
	e2ν = sqrt(e2E*qe*2/source.m)
	        
	sv .= isotropic_velocity(e1ν)

	# assume the new electron and ion are created at the neutral temperature
	T = 300 # target.temperature # 300K = 25°C
	# create new ion and electron
	source.x[source.np+1,:] .= source.x[p,:]
	source.v[source.np+1,:] .= isotropic_velocity(e2ν);
	source.np += 1
	for product in collision.products
		if product == source
			continue
		end
		mp = source.w0/product.w0 + rand()
		νth = thermal_velocity(T, product.m)
		for i=1:round(Integer, mp)
			product.x[product.np+1,:] .= source.x[p,:]
			product.v[product.np+1,:] .= maxwellian_velocity(νth);
			product.np += 1
		end
	end
end

function PIC.perform!(mcc::MonteCarloCollisions, Δt, grid, E)
	ν = zeros(size(grid)) # collision count
	Δh = grid.Δh
	for collision in mcc.collisions
		source, target = collision.source, collision.target
		for p=1:source.np
			i, j, _, _ = PIC.particle_cell(source.x, p, grid.Δh)
			n = PIC.density(target, grid)[i,j]
			if n < 0
				println("Density is negative, skipping")
				continue
			end
			sv = source.v
			tv = (target.q/target.m)*E*Δt
			gv = norm(tv[i,j,:] .- sv[p,:])
			σ  = collision.rate(gv)
			P  = @. 1 - exp(-σ * gv * Δt * n);
			R  = rand()
			if P < R
				continue
			end
			perform!(collision, p, Δt, grid)
			ν[i,j] += 1
		end
	end
	@diag "ν-mcc" PIC.NodeData(ν, grid.origin, [Δh,Δh])
end

function mcc(reactions)
	collisions = []
	for reaction in reactions
		products = []
		source = target = nothing
		@assert length(reaction.reactants) == 2 "Monte Carlo Collisions support only two reacting species: one fluid and one kinetic"
		for (r,_) in reaction.reactants
			if PIC.is_fluid(r) target = r else source = r end
		end
		for (p,c) in reaction.stoichiometry
			if c > 0 push!(products, p) end
		end
		if source ≠ nothing && target ≠ nothing
			if length(products) > 0
				push!(collisions, MCC.IonizationCollision(reaction.rate, source, target, products))
			else
				push!(collisions, MCC.ElasticCollision(reaction.rate, source, target))
			end
		else
			if source == nothing error("Reaction without particle species") end
			if target == nothing error("Reaction without fluid species") end
		end
	end
	MonteCarloCollisions(collisions)
end

Base.show(io :: IO, r :: MCC.ElasticCollision) = print(io, r.rate, ", ", r.source, " + ", r.target, "-->", r.source, " + ", r.target, "\telastic") 
Base.show(io :: IO, r :: MCC.IonizationCollision) = print(io, r.rate, ", ", r.source, " + ", r.target, "-->", r.products, "\tionization") 