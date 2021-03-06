struct ChemicalReactionNetwork
	reactions :: Array{ChemicalReaction,1}
end

function PIC.perform!(network::ChemicalReactionNetwork, E, Δt, config)
	Δn = Dict()
	Δx, Δy = config.grid.Δh
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
		@field "d"*species.name "1/m^2" Δn[species] config.grid
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