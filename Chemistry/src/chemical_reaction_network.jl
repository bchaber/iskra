struct ChemicalReactionNetwork
	reactions :: Array{ChemicalReaction,1}
end

function PIC.perform!(network::ChemicalReactionNetwork, Δt, grid, E)
	Δn = Dict()
	for reaction in network.reactions
		rate = reaction.rate.(E[:,:,1])
		for (species, coeff) in reaction.stoichiometry
			n = species.n
			if haskey(Δn, species) ≠ true
				Δn[species] = zeros(size(n))
			end
			concentration = ones(size(n))
			for (species, order) in reaction.reactants
				concentration .*= species.n.^order
			end
			Δn[species] .+= Δt .* coeff .* rate .* concentration
		end
	end

	for species in keys(Δn)
		species.n .+= Δn[species]
		@diag "Δn"*species.name PIC.NodeData(Δn[species], grid.origin, [1,1]*grid.Δh)
	end
end

function chemical(reactions)
	for reaction in reactions
		for (species, _) in reaction.stoichiometry
			if PIC.is_fluid(species) ≠ true
				error("Chemical reactions are only valid for fluid species: ", species, " is not a fluid.")
			end
		end
	end
	ChemicalReactionNetwork(reactions)
end

Base.show(io :: IO, r :: ChemicalReaction) = print(io, r.rate, ", ", r.reactants, ":", r.stoichiometry) 