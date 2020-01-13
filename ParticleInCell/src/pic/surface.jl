using DataStructures

const BoundaryCell = Tuple{Int64,Int64}
const BoundaryCells = Tuple{Int64,Int64,Int64,Int64}
const TrackedParticle = Tuple{Int64,Int64,Int64,Float64,Float64,Float64}

abstract type Surface end
struct PeriodicSurface <: Surface end
struct AbsorbingSurface <: Surface end
struct ReflectiveSurface <: Surface end

struct SurfaceTracker
  surface :: Dict{BoundaryCells,Surface}
  tracked :: Vector{TrackedParticle}
  Δh :: Float64
end

function Base.get(st::SurfaceTracker,
                  bc::BoundaryCells,
                  dv::Union{Surface,Nothing})
  get(st.surface, bc, dv)
end

function Base.get!(st::SurfaceTracker,
                   bc::BoundaryCells,
                   dv::Union{Surface,Nothing})
  get!(st.surface, bc, dv)
end

function Base.setindex!(st::SurfaceTracker,
                        v0::Surface,
                        bc::BoundaryCells)
  setindex!(st.surface, v0, bc)
end

function Base.in(x::BoundaryCell,  st::SurfaceTracker)
  for (i,j,k,l) in keys(st.surface)
    if (i,j) == x || (k,l) == x
      return true
    end
  end
  return false
end

function Base.in(x::BoundaryCells, st::SurfaceTracker)
  for (i,j,k,l) in keys(st.surface)
    if (i,j,k,l) == x || (k,l,i,j) == x
      return true
    end
  end
  return false
end

create_periodic_surface()   = PeriodicSurface()
create_absorbing_surface()  = AbsorbingSurface()
create_reflective_surface() = ReflectiveSurface()

function reflects(s::Surface)           false end
function reflects(s::ReflectiveSurface) true  end
function absorbs(s::Surface)           true  end
function absorbs(s::PeriodicSurface)   false end
function absorbs(s::ReflectiveSurface) false end
function hit!(s::Surface, part::KineticSpecies, p::Int64,
              Δt::Float64, n̂::Array{Float64,1})
end
function hit!(s::AbsorbingSurface, part::KineticSpecies, p::Int64,
              Δt::Float64, n̂::Array{Float64,1})
end
function hit!(s::ReflectiveSurface, part::KineticSpecies, p::Int64,
              Δt::Float64, n̂::Array{Float64,1})
  part.x[p,:] .-= part.v[p,:]*Δt
  for i=[1,2,3]
    if n̂[i] ≠ 0
      part.v[p,i] *= -1
    end
  end
  part.x[p,:] .+= part.v[p,:]*Δt
end

function build_default_surface!(st::SurfaceTracker, bcs::BitArray{3}, ds::Surface)
  nx, ny = size(bcs)
  for i=1:nx-1
    get!(st, (i,1,i,0),     ds)
    get!(st, (i,ny-1,i,ny), ds)
  end

  for j=1:ny-1
    get!(st, (1,j,0,j),     ds)
    get!(st, (nx-1,j,nx,j), ds)
  end
end

function build_surface_lookup!(st::SurfaceTracker, bcs::BitArray{3}, ss::Surface)
  nx, ny = size(bcs)
  for i=1:nx-1
    for j=1:ny-1
      a = bcs[i,  j  ,1]
      b = bcs[i+1,j  ,1]
      c = bcs[i+1,j+1,1]
      d = bcs[i  ,j+1,1]
      if a st[(i,j, (a == b) ? i : i-1, j-1)] = ss end
      if b st[(i,j, i+1, (b == c) ? j : j-1)] = ss end
      if c st[(i,j, (c == d) ? i : i+1, j+1)] = ss end
      if d st[(i,j, i-1, (d == a) ? j : j+1)] = ss end
    end
  end
end

function create_surface_tracker(grid::UniformGrid{XY2D})
  Δx, _, _ = grid.Δh
  nx, ny = size(grid)
  bcs = zeros(Int64, nx, ny, 1)

  tracked = TrackedParticle[]
  surface = Dict{BoundaryCells, Surface}()
  ds = AbsorbingSurface()
  st = SurfaceTracker(surface, tracked, Δx)

  build_default_surface!(st, bcs .== 0, ds)
  return st
end

function track_surface!(st::SurfaceTracker, bcs::BitArray{3}, ss::Surface)
  build_surface_lookup!(st, bcs, ss)
end

function create_surface_tracker(bcs::Array{Int8, 3}, ss::Array{<:Surface,1}, Δh)
  tracked = TrackedParticle[]
  surface = Dict{BoundaryCells, Surface}()
  st = SurfaceTracker(surface, tracked, Δh)
  ds = AbsorbingSurface()

  build_default_surface!(st, bcs .== 0, ds)
  for i=1:length(ss)
    build_surface_lookup!(st, bcs .== i, ss[i])
  end
  return st
end

function track!(::Nothing, part::KineticSpecies, Δt) end
function track!(st::SurfaceTracker, part::KineticSpecies, Δt)
  Δh = st.Δh
  np = part.np
  px = view(part.x, 1:np, :)
  pv = view(part.v, 1:np, :)

  deleteat!(st.tracked, 1:length(st.tracked))
  for p=1:part.np
    i, j, hx, hy = particle_cell(px, p, Δh)
    if (i,j) ∈ st
      push!(st.tracked, (p,i,j,hx,hy,Δt))
    end
  end
end

function check_particle(i, j, hx, hy, Δt, Δh, vx, vy)
  dx = (vx > 0) ? Δh*(1-hx) : Δh*hx
  dy = (vy > 0) ? Δh*(1-hy) : Δh*hy

  Δtx, Δty = dx/abs(vx), dy/abs(vy)
  if Δt < Δtx && Δt < Δty # particle stayed in the same cell
    return i, j, hx, hy, Δt
  end

  if Δtx < Δty
    return (vx > 0) ?
      (i+1, j, 0., hy + vy*Δtx/Δh, Δt-Δtx) :
      (i-1, j, 1., hy + vy*Δtx/Δh, Δt-Δtx)
  else
    return (vy > 0) ?
      (i, j+1, hx + vx*Δty/Δh, 0., Δt-Δty) :
      (i, j-1, hx + vx*Δty/Δh, 1., Δt-Δty)
  end
end
    
function check!(::Nothing, part::KineticSpecies, Δt) end
function check!(st::SurfaceTracker, part::KineticSpecies, Δt)
  Δh = st.Δh
  px = view(part.x, 1:part.np, :)
  pv = view(part.v, 1:part.np, :)
  absorbed = SortedSet{Int64}(Base.Order.Reverse)

  vmax = Δh./Δt
  pvmax = part.np > 0 ? maximum(abs.(pv); dims=1) : 0.0
  if any(pvmax .> vmax)
    println("ERROR: $part particle is too fast: $pvmax > $vmax")
  end

  while length(st.tracked) > 0
    p, i , j , hx , hy , Δt  = popfirst!(st.tracked)
    vx, vy = pv[p,1], pv[p,2]
    
    i′, j′, hx′, hy′, Δt′ = check_particle(i, j, hx, hy, Δt, Δh, vx, vy)
    if i == i′ && j == j′
      continue
    end
    n̂ = [i′-i, j′-j, 0.]
    surface = get(st, (i, j, i′,j′), nothing)
    if surface ≠ nothing
      if absorbs(surface)
        push!(absorbed, p)
      end

      if reflects(surface)
        γ = (Δt-Δt′)/Δh
        push!(st.tracked, (p, i, j, hx+γ*vx, hy+γ*vy, Δt′))
      end

      hit!(surface, part, p, Δt′, n̂)
    else
      push!(st.tracked, (p, i′, j′, hx′, hy′, Δt′))
    end
  end

  for p in absorbed
    remove!(part, p)
  end
end

function wrap!(part::KineticSpecies, grid::UniformGrid{XY2D})
  nx, ny, nz = grid.nx-1, grid.ny-1, grid.nz-1
  Δx, Δy, Δz = grid.Δh
  Lx, Ly, Lz = [nx, ny, nz] .* [Δx, Δy, Δz]
  origin = grid.origin
  bb = [origin origin] .+ [0. Lx; 0. Ly; 0 Lz]
  px = view(part.x, 1:part.np, :)
  for p=1:part.np
    if px[p,1] < bb[1,1] px[p,1] += Lx end
    if px[p,1] > bb[1,2] px[p,1] -= Lx end
    if px[p,2] < bb[2,1] px[p,2] += Ly end
    if px[p,2] > bb[2,2] px[p,2] -= Ly end
  end
end