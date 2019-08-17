module Pusher
export push_particles!
struct BorisPusher end

function create_boris_pusher()
  return BorisPusher()
end

function push_particles!(::BorisPusher, part, E, Δt)
  np, id = part.np, part.id
  x, v, q, m = part.x, part.v, part.q, part.m
  F = q*E # Lorentz force, F=qE
  a = F/m # acceleration,  F=ma

  v[1:np,:] .= @views v[1:np,:] .+ Δt * a[1:np,:]
  x[1:np,:] .= @views x[1:np,:] .+ Δt * v[1:np,:]
end
end
