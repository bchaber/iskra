
function check(pt::TrackedParticle, pv::AbstractArray{Float64,2}, Δh)
  Δt, p, i, j, hx, hy = pt
  vx, vy = pv[p,1], pv[p,2]
  dx = (vx > 0) ? Δh*(1-hx) : Δh*hx
  dy = (vy > 0) ? Δh*(1-hy) : Δh*hy

  Δtx, Δty = dx/abs(vx), dy/abs(vy)
  if Δt < Δtx && Δt < Δty # particle stayed in the same cell
    return pt
  end

  if Δtx < Δty
    return (vx > 0) ?
      (Δt-Δtx, p, i+1, j, 0., hy + vy*Δtx/Δh) :
      (Δt-Δtx, p, i-1, j, 1., hy + vy*Δtx/Δh)
  else
    return (vy > 0) ?
      (Δt-Δty, p, i, j+1, hx + vx*Δty/Δh, 0.) :
      (Δt-Δty, p, i, j-1, hx + vx*Δty/Δh, 1.)
  end
end
    
function check!(::Nothing, part::KineticSpecies, Δt) end
function check!(st::SurfaceTracker, part::KineticSpecies, Δt)
  # Check if there are particles moving too fast
  pv = view(part.v, 1:part.np, :)
  vmax = st.Δh./Δt
  pvmax = part.np > 0 ? maximum(abs.(pv); dims=1) : 0.0
  if any(pvmax .> vmax)
    println("ERROR: $part particle is too fast: $pvmax > $vmax")
  end
  empty!(st.absorbed)
  while length(st.tracked) > 0
    pt  = popfirst!(st.tracked)
    pt′ = check(pt, pv, st.Δh)
    _, p, i,  j  = pt
    _, _, i′, j′ = pt′
    if i == i′ && j == j′
      continue
    end
    surface = get(st, (i, j, i′,j′), nothing)
    if surface ≠ nothing
      hit!(surface, part, st, pt, pt′)
    else
      track!(st, pt′)
    end
  end
  # Remove absorbed particles
  for p in st.absorbed
    remove!(part, p)
  end
  # Emit new particles
  # ...
end