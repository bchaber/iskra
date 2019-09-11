mutable struct FluidSpecies
  name::String
	μ :: Float64 # mobility of the species
	q :: Float64 # charge of a single particle
	m :: Float64 # mass of a single particle
	n :: Array{Float64,2} # density of the species
end

density(species :: FluidSpecies, grid) = species.n