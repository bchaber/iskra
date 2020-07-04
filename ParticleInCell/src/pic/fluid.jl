mutable struct FluidSpecies{D}
  name::String
	μ :: Float64 # mobility of the species
	q :: Float64 # charge of a single particle
	m :: Float64 # mass of a single particle
	n :: Array{Float64,D} # density of the species
	T :: Float64 # temperature
end
Base.show(io::IO, sp::FluidSpecies) = print(io, sp.name)
is_fluid(species :: FluidSpecies) = true
density(species :: FluidSpecies, _) = species.n