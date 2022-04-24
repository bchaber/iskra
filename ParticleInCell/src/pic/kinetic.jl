mutable struct KineticSpecies{D, V}
  x :: Vector{SVector{D, Float64}} # position
  v :: Vector{SVector{V, Float64}} # velocity
  m :: Float64 # mass
  q :: Float64 # charge
  np :: Int64  # current number of particles
  w0 :: Float64 # default weight of new particle
  wg :: Vector{Float64} # number of particles per macroparticle
  id :: Vector{UInt32} # identifier
  name :: Symbol
end

Base.show(io::IO, sp::KineticSpecies) = print(io, sp.name)
particle_uuids(N::Int64) = collect(UInt32(1):UInt32(N))
KineticSpecies{D,V}(name::Symbol, N::Int64, q=0.0, m=0.0, weight=1.0) where {D, V} =
  KineticSpecies{D,V}(
    zeros(SVector{D, Float64}, N),
    zeros(SVector{V, Float64}, N),
    m, q, 0, weight, weight*ones(N), particle_uuids(N), name)
function remove!(sp::KineticSpecies, i::Int64)
  np = sp.np
  sp.x[i] = sp.x[np]
  sp.v[i] = sp.v[np]
  sp.wg[i], sp.wg[np] = sp.wg[np], sp.w0
  sp.id[i], sp.id[np] = sp.id[np], sp.id[i]
  sp.np = np - 1
end

function add!(src, dst)
  if src.np > 0
    rn = (1 + dst.np):(dst.np + src.np)
    dst.x[rn] = src.x
    dst.v[rn] = src.v
    # dst has already correct ID numbers
    dst.np += src.np
  end
end

function remove_particles!(part, Δh, matches)
  p = 1
  while p <= part.np
    i, j, _, _ = particle_cell(part, p, Δh)
    if matches(i,j)
      remove!(part, p)
      continue
    end
    p = p + 1
  end
end

is_fluid(species :: KineticSpecies) = false