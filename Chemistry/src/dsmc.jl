using LinearAlgebra
import ParticleInCell

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
	collision_count :: Array{Float64,2}
	collision_source_candidates :: Array{Array{UInt32,1},2}
	collision_target_candidates :: Array{Array{UInt32,1},2}
end

DirectSimulationMonteCarlo(collisions) = DirectSimulationMonteCarlo(collisions, zeros(0,0), zeros(0,0), zeros(0,0))

function cache!(candidates, species, nx, ny, Δh)
	for p=1:species.np
		i, j, ~, ~ = ParticleInCell.particle_cell(species.x, p, Δh)
		push!(candidates[i, j], p)
	end
end

function perform!(collision::DSMC.ElasticCollision, s, t)
	source, target = collision.source, collision.target
	
	mr1 = source.m/(source.m + target.m)
	mr2 = target.m/(source.m + target.m)

	g = (source.v[s,:] - target.v[t,:])
	vc_cm = mr1 * source.v[s,:] .+ mr2 * target.v[t,:]
	vss_inv= 1.
	print("*")
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
	Δh = config.grid.Δh
	ν = zeros(nx, ny) # collision count
	dsmc.collision_source_candidates = [Vector{UInt32}() for i=1:nx, j=1:ny]
	dsmc.collision_target_candidates = [Vector{UInt32}() for i=1:nx, j=1:ny]
	for collision in dsmc.collisions
		source, target = collision.source, collision.target
		cache!(dsmc.collision_source_candidates, source, nx, ny, Δh) # assign particles to cells
		cache!(dsmc.collision_target_candidates, target, nx, ny, Δh) # assign particles to cells
		Pab = target.w0/source.w0
		Pba = source.w0/target.w0

		if Pab > 1 Pab = 1 end
		if Pba > 1 Pba = 1 end
		
		σmax = maximum(collision.rate)
		N  = source.np * target.np * Δt/Δh^2 * σmax
		N /= Pab + 1
		N = floor(Int64, N)
		for i=1:nx
			for j=1:ny
				for n=1:N
					if length(dsmc.collision_source_candidates[i,j]) == 0 || 
					   length(dsmc.collision_target_candidates[i,j]) == 0
						continue
					end

				    sR = rand(dsmc.collision_source_candidates[i,j])
					tR = rand(dsmc.collision_target_candidates[i,j])
					while source ≠ target && sR == tR
						tR = rand(dsmc.collision_target_candidates[i,j])
					end

					cr = norm(source.v[sR,:] - target.v[tR,:])
					σcr = collision.rate(cr)

				    P = σcr / σmax
				    R = rand()
				    
				    if P > R
						perform!(collision, sR, tR)
						ν[i,j] += 1
					end
				end
			end
		end
	end
	@diag "ν" PIC.NodeData(ν, config.grid.origin, [Δh,Δh])
end

function dsmc(reactions)
	reminder = 0
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