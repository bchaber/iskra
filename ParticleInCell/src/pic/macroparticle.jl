mutable struct Species
  x :: AbstractArray{Float64,2} # position
  v :: AbstractArray{Float64,2} # velocity
  m :: Float64 # mass
  q :: Float64 # charge
  name::String
  np2c :: Int64 # number of particles per macroparticle
  np::Int64     # current number of particles
end

Species(name::String, N::Int64) = Species(zeros(N,2), zeros(N,2), 0, 1, name, 1, 0)

function remove!(sp, i::Int64)
  np = sp.np
  sp.x[i,:] .= sp.x[np,:]
  sp.v[i,:] .= sp.v[np,:]
  sp.np = np - 1
end

function add!(src, dst)
  rn = (1 + dst.np):(dst.np + src.np)
  dst.x[rn,:] .= src.x
  dst.v[rn,:] .= src.v
  dst.np += src.np
end
