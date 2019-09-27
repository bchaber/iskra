module Chemistry

export @reactions
export CrossSection
export ChemicalReaction
export ChemicalReactionNetwork
export MonteCarloCollisions
export DirectSimulationMonteCarlo
export ElasticCollision, IonizationCollision
export DSMCElasticCollision, DSMCIonizationCollision
export chemical, mcc, dsmc

is_fluid(x) = error("Implement is_fluid for ", typeof(x))

include("reactions.jl")
include("cross_section.jl")
include("chemical_reaction_network.jl")
include("mcc.jl")
include("dsmc.jl")

end