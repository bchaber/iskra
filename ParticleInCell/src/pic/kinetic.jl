mutable struct KineticSpecies{D, V}
  x :: AbstractArray{Float64,2} # position
  v :: AbstractArray{Float64,2} # velocity
  n :: AbstractArray{Float64,2} # density
  m :: Float64 # mass
  q :: Float64 # charge
  np :: Int64     # current number of particles
  w0 :: Float64   # default weight of new particle
  wg :: AbstractArray{Float64,1}  # number of particles per macroparticle
  id :: AbstractArray{UInt32,1} # identifier
  name::String
end

Base.show(io::IO, sp::KineticSpecies) = print(io, sp.name)
particle_uuids(N::Int64) = collect(UInt32(1):UInt32(N))
KineticSpecies{D,V}(name::String, N::Int64) where {D, V} =
  KineticSpecies{D,V}(zeros(N, D), zeros(N, V), zeros(0,0),
    .0, .0, 0, 1., ones(N), particle_uuids(N), name)

function remove!(sp::KineticSpecies, i::Int64)
  np = sp.np
  sp.x[i,:] .= sp.x[np,:]
  sp.v[i,:] .= sp.v[np,:]
  sp.wg[i], sp.wg[np] = sp.wg[np], sp.w0
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

function remove_particles!(part, Δh, matches)
  p = 1
  px = view(part.x, 1:part.np, :)
  while p <= part.np
    i, j, _, _ = particle_cell(px, p, Δh)
    if matches(i,j)
      remove!(part, p)
      continue
    end
    p = p + 1
  end
end

is_fluid(species :: KineticSpecies) = false
density(species :: KineticSpecies, grid) = particle_to_grid(species, grid, (p) -> species.wg[p])