module Chemistry

export @reactions
export CrossSection
export Reaction
export react

include("reactions.jl")
include("cross_section.jl")

Base.show(io :: IO, r :: Reaction) = print(io, r.rate, ", ", r.reactants, ":", r.stoichiometry) 
end