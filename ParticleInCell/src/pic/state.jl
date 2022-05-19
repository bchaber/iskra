struct ParticleInCellState{D, V, F, S}
    particles :: NTuple{S, KineticSpecies{D, V}}
    potential :: Array{F, D}
    electric :: Array{SVector{V, F}, D}
    magnetic :: Array{SVector{V, F}, D}
    density :: Array{F, D}
    current :: Array{F, D}

    timestep :: Float64
    cellsize :: Float64
end