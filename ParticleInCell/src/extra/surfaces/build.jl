using DataStructures

const BoundaryCell{D} = NTuple{D, Int64}
const BoundaryCells{D} = Tuple{BoundaryCell{D}, BoundaryCell{D}}
const TrackedParticle{D} = Tuple{Float64, Int64, NTuple{D,Int64}, NTuple{D, Float64}}
const AbsorbedParticle = Int64

abstract type Surface end
struct PeriodicSurface <: Surface end
struct AbsorbingSurface <: Surface end
struct ReflectiveSurface <: Surface end

struct SurfaceTracker{D}
  surface :: Dict{BoundaryCells{D}, Surface}
  tracked :: Vector{TrackedParticle{D}}
  absorbed :: SortedSet{AbsorbedParticle}
  Δh :: Float64
end

export create_periodic_surface
export create_absorbing_surface
export create_reflective_surface
export track_surface!
create_periodic_surface()   = PeriodicSurface()
create_absorbing_surface()  = AbsorbingSurface()
create_reflective_surface() = ReflectiveSurface()

function build_default_surface!(st::SurfaceTracker{1}, bcs::BitArray{1}, ds::Surface)
  nx, = size(bcs)
  get!(st,    ((1,),  (0,)), ds)
  get!(st, ((nx-1,), (nx,)), ds)
end

function build_default_surface!(st::SurfaceTracker{2}, bcs::BitArray{2}, ds::Surface)
  nx, ny = size(bcs)
  for i=1:nx-1
    get!(st, ((i,    1), (i,  0)), ds)
    get!(st, ((i, ny-1), (i, ny)), ds)
  end

  for j=1:ny-1
    get!(st, ((   1, j), ( 0, j)), ds)
    get!(st, ((nx-1, j), (nx, j)), ds)
  end
end

function build_surface_lookup!(st::SurfaceTracker{2}, bcs::BitArray{2}, ss::Surface)
  nx, ny = size(bcs)
  for i=1:nx-1
    for j=1:ny-1
      a = bcs[i,  j  ,1]
      b = bcs[i+1,j  ,1]
      c = bcs[i+1,j+1,1]
      d = bcs[i  ,j+1,1]
      if (a == b == true) st[((i,j), (i, j-1))] = ss end
      if (b == c == true) st[((i,j), (i+1, j))] = ss end
      if (c == d == true) st[((i,j), (i, j+1))] = ss end
      if (d == a == true) st[((i,j), (i-1, j))] = ss end
    end
  end
end

function create_surface_tracker(bcs::Array{Int8,1}, ss::Vector{<:Surface}, Δh,
        ds = AbsorbingSurface())
  tracked = Vector{TrackedParticle{1}}()
  surface = Dict{BoundaryCells{1}, Surface}()
  absorbed = SortedSet{AbsorbedParticle}(Base.Order.Reverse)
  st = SurfaceTracker{1}(surface, tracked, absorbed, Δh)

  build_default_surface!(st, bcs .== 0, ds)
  for i=1:length(ss)
    build_surface_lookup!(st, bcs .== i, ss[i])
  end
  return st
end

function create_surface_tracker(bcs::Array{Int8,2}, ss::Vector{<:Surface}, Δh,
        ds = AbsorbingSurface())
  tracked = Vector{TrackedParticle{2}}()
  surface = Dict{BoundaryCells{2}, Surface}()
  absorbed = SortedSet{AbsorbedParticle}(Base.Order.Reverse)
  st = SurfaceTracker{2}(surface, tracked, absorbed, Δh)

  build_default_surface!(st, bcs .== 0, ds)
  for i=1:length(ss)
    build_surface_lookup!(st, bcs .== i, ss[i])
  end
  return st
end

function create_surface_tracker(grid::CartesianGrid{2}, ds=AbsorbingSurface())
  Δx, ~  = grid.Δh
  nx, ny = size(grid)
  bcs = zeros(Int8, nx, ny)
  return create_surface_tracker(bcs, Surface[], Δx, ds)
end

function create_surface_tracker(grid::CartesianGrid{1}, ds=AbsorbingSurface())
  Δx, = grid.Δh
  nx, = size(grid)
  bcs = zeros(Int8, nx)
  return create_surface_tracker(bcs, Surface[], Δx, ds)
end

function track_surface!(st::SurfaceTracker, bcs::BitArray{2}, ss::Surface)
  build_surface_lookup!(st, bcs, ss)
end
