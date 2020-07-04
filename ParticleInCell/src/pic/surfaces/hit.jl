function normalvector(pt :: TrackedParticle, pt′ :: TrackedParticle)
  _, _, i,  j  = pt
  _, _, i′, j′ = pt′
  return [i′-i, j′-j, 0.]
end

function absorbed!(st::SurfaceTracker, pt′::TrackedParticle)
  _, p, _ = pt′
  push!(st.absorbed, p)
end

function scattered!(st::SurfaceTracker, pt′::TrackedParticle)
  Δt, p, i, j, hx, hy = pt′
  i′, j′, hx′, hy′ = i, j, hx, hy
  if (hx == 0.)  i′, hx′ = i-1, 1. end
  if (hx == 1.)  i′, hx′ = i+1, 0. end
  if (hy == 0.)  j′, hy′ = j-1, 1. end
  if (hy == 1.)  j′, hy′ = j+1, 0. end
  track!(st, (Δt, p, i′, j′, hx′, hy′))
end

function emit!(st::SurfaceTracker, pt′::TrackedParticle,
               v′::AbstractVector{Float64})
end

function hit!(::Surface,
              ::KineticSpecies,
              ::SurfaceTracker,
              pt ::TrackedParticle,
              pt′::TrackedParticle)
end
function hit!(::AbsorbingSurface,
              ::KineticSpecies,
              st ::SurfaceTracker,
              pt ::TrackedParticle,
              pt′::TrackedParticle)
  absorbed!(st, pt′)
end
function hit!(s::ReflectiveSurface,
              part::KineticSpecies{D,V},
              st ::SurfaceTracker,
              pt ::TrackedParticle,
              pt′::TrackedParticle) where {D,V}
  Δt′, p, i′, j′ = pt′
  n̂  = normalvector(pt, pt′)
  px = view(part.x, p, 1:D)
  pv = view(part.v, p, 1:D)

  px .-= pv*Δt′
  for i=1:3
    if n̂[i] ≠ 0
      pv[i] *= -1
    end
  end
  px .+= pv*Δt′
  scattered!(st, pt′)
end