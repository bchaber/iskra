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

struct DirectSimulationMonteCarlo
	collisions
	cells
	reminder
end

function PIC.perform!(dsmc::DirectSimulationMonteCarlo, E, Δt, config)
	ν = zeros(size(config.grid)) # collision count
	Δh = config.grid.Δh
	for collision in dsmc.collisions
		source, target = collision.source, collision.target
		cache(source, target, dsmc) # assign particles to cells
	end
	@diag "ν-mcc" PIC.NodeData(ν, config.grid.origin, [Δh,Δh])
end

function dsmc(reactions)
	cells = Dict{Tuple{Int64,Int64}, Array{Int64,1}}()
	reminder = 0
	collisions = []
	for reaction in reactions
		products = []
		@assert length(reaction.reactants) == 2 "Direct Simulation Monte Carlo support only two reacting, kinetic species"
		source = reaction.reactants[1]
		target = reaction.reactants[2]
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
	DirectSimulationMonteCarlo(collisions, cells, reminder)
end