
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

function Base.in(x::BoundaryCell,  st::SurfaceTracker{D}) where {D} 
  for (ij, kl) in keys(st.surface)
    if ij == x || kl == x
      return true
    end
  end
  return false
end

function Base.in(x::BoundaryCells{D}, st::SurfaceTracker{D}) where {D}
  for (ij, kl) in keys(st.surface)
    if (ij..., kl...) == x || (kl..., ij...) == x
      return true
    end
  end
  return false
end

function track!(st::SurfaceTracker, pt::TrackedParticle)
  push!(st.tracked, pt)
end
function track!(  ::Nothing,        part::KineticSpecies, Δt) end
function track!(st::SurfaceTracker, part::KineticSpecies, Δt)
  px = view(part.x, 1:part.np, :)

  empty!(st.tracked)
  for p=1:part.np
    i, j, hx, hy = particle_cell(px, p, st.Δh)
    if (i, j) ∈ st
      track!(st, (Δt, p, (i, j), (hx, hy)))
    end
  end
end
