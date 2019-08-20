mutable struct Species
  x :: AbstractArray{Float64,2} # position
  v :: AbstractArray{Float64,2} # velocity
  m :: Float64 # mass
  q :: Float64 # charge
  name::String
  np2c :: Int64 # number of particles per macroparticle
  np :: Int64     # current number of particles
  id :: AbstractArray{UInt32,1}   # identifier
end

particle_uuids(N::Int64) = collect(1:UInt32(N))

Species(x::AbstractArray{Float64,2},
        v::AbstractArray{Float64,2},
        m::Float64, q::Float64, name::String,
        np2c::Int64, np::Int64) = Species(x, v, m, q, name, np2c, np, particle_uuids(np))
Species(name::String, N::Int64) = Species(zeros(N,3), zeros(N,3), 0, 1, name, 1, 0, particle_uuids(N))

function remove!(sp, i::Int64)
  np = sp.np
  sp.x[i,:] .= sp.x[np,:]
  sp.v[i,:] .= sp.v[np,:]
  sp.id[i], sp.id[np] = sp.id[np], sp.id[i]
  sp.np = np - 1
end

function add!(src, dst)
  if src.np > 0
    rn = (1 + dst.np):(dst.np + src.np)
    dst.x[rn,:] .= src.x
    dst.v[rn,:] .= src.v
    # dst has already correct ID numbers
    dst.np += src.np
  end
end
