module Chemistry

export @reactions
export CrossSection
export ChemicalReaction
export ChemicalReactionNetwork
export chemical, mcc, dsmc

is_fluid(x) = error("Implement is_fluid for ", typeof(x))

import LinearAlgebra: norm, dot
import Diagnostics: @field
import ParticleInCell
PIC = ParticleInCell

include("reactions.jl")
include("cross_section.jl")
include("chemical_reaction_network.jl")
include("see.jl")
include("mcc.jl")
include("dsmc.jl")
include("lxcat.jl")
end