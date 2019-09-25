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

struct MonteCarloCollisions
	collisions
end

function mcc(reactions)
	collisions = []
	for reaction in reactions
		products = []
		source = target = nothing
		@assert length(reaction.reactants) == 2 "Monte Carlo Collisions support only two reacting species: one fluid and one kinetic"
		for (r,_) in reaction.reactants
			if is_fluid(r) target = r else source = r end
		end
		for (p,c) in reaction.stoichiometry
			if c > 0 push!(products, p) end
		end
		if source ≠ nothing && target ≠ nothing
			if length(products) > 0
				push!(collisions, IonizationCollision(reaction.rate, source, target, products))
			else
				push!(collisions, ElasticCollision(reaction.rate, source, target))
			end
		else
			if source == nothing error("Reaction without particle species") end
			if target == nothing error("Reaction without fluid species") end
		end
	end
	MonteCarloCollisions(collisions)
end

Base.show(io :: IO, r :: ElasticCollision) = print(io, r.rate, ", ", r.source, " + ", r.target, "-->", r.source, " + ", r.target, "\telastic") 
Base.show(io :: IO, r :: IonizationCollision) = print(io, r.rate, ", ", r.source, " + ", r.target, "-->", r.products, "\tionization") 