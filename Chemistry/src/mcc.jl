struct Collision
	rate
	source # particles
	target # fluid
	products # particles
end

struct MonteCarloCollisions
	collisions :: Array{Collision,1}
end

function mcc(reactions)
	collisions = Collision[]
	for reaction in reactions
		products = []
		source = target = nothing
		for (r,_) in reaction.reactants
			if is_fluid(r) target = r else source = r end
		end
		for (p,c) in reaction.stoichiometry
			if c > 0 push!(products, p) end
		end
		if source ≠ nothing && target ≠ nothing
			push!(collisions, Collision(reaction.rate, source, target, products))
		else
			if source == nothing error("Reaction without particle species") end
			if target == nothing error("Reaction without fluid species") end
		end
	end
	MonteCarloCollisions(collisions)
end

Base.show(io :: IO, c :: Collision) = print(io, c.source, "->", c.target, ":", c.products) 