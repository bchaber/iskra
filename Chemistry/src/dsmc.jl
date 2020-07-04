module DSMC
	struct ElasticCollision
		rate
		source
		target
	end

	struct IonizationCollision
		rate
		source
		target
		products
	end
end

mutable struct DirectSimulationMonteCarlo
	collisions
	collisions_remaining :: Array{Float64,2}
	collision_source_candidates :: Array{Array{UInt32,1},2}
	collision_target_candidates :: Array{Array{UInt32,1},2}
end

DirectSimulationMonteCarlo(collisions) = DirectSimulationMonteCarlo(collisions, zeros(0,0), zeros(0,0), zeros(0,0))

function cache!(candidates, species, nx, ny, Δh)
	for p=1:species.np
		i, j, ~, ~ = PIC.particle_cell(species.x, p, Δh)
		push!(candidates[i, j], p)
	end
end

function perform!(collision::DSMC.ElasticCollision, s, t)
	source, target = collision.source, collision.target
	
	mr1 = source.m/(source.m + target.m)
	mr2 = target.m/(source.m + target.m)

	g = source.v[s,:] .- target.v[t,:]
	vc_cm = mr1 * source.v[s,:] .+ mr2 * target.v[t,:]
	vss_inv= 1.
	#print("*")
	if abs(vss_inv - 1.0) < 1e-4
		B = 2rand() - 1.0
		A = sqrt(1 - B^2)
		C = 2π*rand()
		vr_cp = norm(g) * [B, A*cos(C), A*sin(C)]
	else
		B = 2rand()^(vss_inv) - 1
		A = sqrt(1 - B^2)
		C = 2π*rand()
		D = sqrt(dot(g, g))

		if D > 1e-6
			vr_cp = B * g
			vr_cp[1] += (A*D)*sin(C)
			vr_cp[2] += (A/D)*(norm(g)*cos(C))*g[3] - (A/D)*(g[1]*sin(C))*g[2]
			vr_cp[3] -= (A/D)*(norm(g)*cos(C))*g[2] - (A/D)*(g[1]*sin(C))*g[3]
		else
			vr_cp = g[1]*[B, A*cos(C), A*sin(C)]
		end
	end

	if source.wg[s] == target.wg[t]
		source.v[s,:] = vc_cm .+ mr2 * vr_cp
		target.v[t,:] = vc_cm .- mr1 * vr_cp
	else
		Pab = target.wg[t]/source.wg[s]
		Pba = source.wg[s]/target.wg[t]
		R = rand()

		if Pab > R
			source.v[s,:] = vc_cm .+ mr2 * vr_cp
		end

		if Pba > R
			target.v[s,:] = vc_cm .- mr2 * vr_cp
		end
	end
end

function PIC.perform!(dsmc::DirectSimulationMonteCarlo, E, Δt, config)
	nx, ny = size(config.grid)
	Δx, Δy = config.grid.Δh
	ν = zeros(nx, ny) # collision count
	if length(dsmc.collisions_remaining) == 0
		dsmc.collisions_remaining = zeros(nx, ny)
	end
	dsmc.collision_source_candidates = [Vector{UInt32}() for i=1:nx, j=1:ny]
	dsmc.collision_target_candidates = [Vector{UInt32}() for i=1:nx, j=1:ny]
	for collision in dsmc.collisions
		source, target = collision.source, collision.target
		cache!(dsmc.collision_source_candidates, source, nx, ny, config.grid.Δh) # assign particles to cells
		cache!(dsmc.collision_target_candidates, target, nx, ny, config.grid.Δh) # assign particles to cells
		Wa, Wb = source.w0, target.w0
		if Wa > Wb
			Pab, Pba = Wb/Wa, 1.
		else
			Pab, Pba = 1., 	Wa/Wb
		end
		
		σgmax = maximum(collision.rate) * argmax(collision.rate)
		for i=1:nx
			for j=1:ny
				Na = length(dsmc.collision_source_candidates[i,j])
				Nb = length(dsmc.collision_target_candidates[i,j])
				if Na < 2 || Nb < 2
					continue
				end
				na  = Na * Wa / (Δx*Δy)
				Nc  = na * Nb * Δt * σgmax
				Nc /= Pab + (Wb/Wa) * Pba
				if source ≠ target
					Nc *= 2
				end

				Nc += dsmc.collisions_remaining[i,j]
				dsmc.collisions_remaining[i,j] = Nc - floor(Int64, Nc)
				
				for ~=1:floor(Int64, Nc)
				    sR = rand(dsmc.collision_source_candidates[i,j])
					tR = rand(dsmc.collision_target_candidates[i,j])
					while source ≠ target && sR == tR
						tR = rand(dsmc.collision_target_candidates[i,j])
					end

					g = norm(source.v[sR,:] - target.v[tR,:])
					σg = collision.rate(g) * g

				    P = σg / σgmax
				    R = rand()
				    
				    if P < R
				    	continue
				    end
					perform!(collision, sR, tR)
					ν[i,j] += 1
				end
			end
		end
	end
	@field "nuDSMC" "1/m^2" ν config.grid
end

function dsmc(reactions)
	reminder = 0.0
	collisions = []
	for reaction in reactions
		products = []
		@assert length(reaction.reactants) == 2 "Direct Simulation Monte Carlo support only two reacting, kinetic species"
		source, _ = reaction.reactants[1]
		target, _ = reaction.reactants[2]
		for (p,c) in reaction.stoichiometry
			if c > 0 push!(products, p) end
		end
		if source ≠ nothing && target ≠ nothing
			if length(products) > 0
				push!(collisions, DSMC.IonizationCollision(reaction.rate, source, target, products))
			else
				push!(collisions, DSMC.ElasticCollision(reaction.rate, source, target))
			end
		else
			if source == nothing error("Reaction without source species") end
			if target == nothing error("Reaction without target species") end
		end
	end
	DirectSimulationMonteCarlo(collisions)
end