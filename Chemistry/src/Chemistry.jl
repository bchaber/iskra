module Chemistry

export @reactions
export CrossSection
export Reaction
export react!

include("reactions.jl")
include("cross_section.jl")

function react!(x)
	return 0
end

Base.show(io :: IO, r :: Reaction) = print(r.rate, ", ", r.reactants, ":", r.stoichiometry) 
end