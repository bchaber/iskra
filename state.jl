struct ParticleInCellState{D, V, F, S} where {D <: Number, V <: Number, S <: Number, F}
    particles :: NTuple{KineticSpecies{D, V}, S}
    potential :: Array{F, D}
    electric :: Array{SVector{V, F}, D}
    magnetic :: Array{SVector{V, F}, D}
    density :: Array{F, D}
    current :: Array{F, D}

    timestep :: Float64
    cellsize :: Float64
end

struct TwoStreamInstability

end