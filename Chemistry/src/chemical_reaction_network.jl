struct ChemicalReactionNetwork
	reactions :: Array{ChemicalReaction,1}
end

function chemical(reactions)
	for reaction in reactions
		for (species, _) in reaction.stoichiometry
			if is_fluid(species) ≠ true
				error("Chemical reactions are only valid for fluid species: ", species, " is not a fluid.")
			end
		end
	end
	ChemicalReactionNetwork(reactions)
end

Base.show(io :: IO, r :: ChemicalReaction) = print(io, r.rate, ", ", r.reactants, ":", r.stoichiometry) 