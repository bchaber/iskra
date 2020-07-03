export create_boris_pusher

struct BorisPusher end
function create_boris_pusher()
  return BorisPusher()
end

function push_particles!(::BorisPusher,
	part::KineticSpecies{D, V}, E, Δt) where {D, V, T}
  np, id = part.np, part.id
  x, v, q, m = part.x, part.v, part.q, part.m
  F = q*E # Lorentz force, F=qE
  a = F/m # acceleration,  F=ma

  v[1:np,1:V] .= @views v[1:np,1:V] .+ Δt * a[1:np,1:V]
  x[1:np,1:D] .= @views x[1:np,1:D] .+ Δt * v[1:np,1:D]
end