module Pusher
export push_particles!
import Diagnostics
struct BorisPusher end

function create_boris_pusher()
  return BorisPusher()
end

struct ParticleVectorData <: Diagnostics.DiagnosticData
  x :: Array{Float64,1}
  y :: Array{Float64,1}
  u :: Array{Float64,1}
  v :: Array{Float64,1}
  w :: Array{Float64,1}
 id :: Array{UInt64,1}
 it :: Integer
end

import PlotVTK: pvd_add_timestep, scatter
Diagnostics.save_diagnostic(dname::String, d::ParticleVectorData, cname::String, c::Any, it::Integer) =
  pvd_add_timestep(c, scatter(d.x, d.y, dname, dname => (d.u, d.v, d.w), "uuid" => (d.id,), it=it, save=false), it)

function push_particles!(::BorisPusher, part, E, Δt)
  np, id = part.np, part.id
  x, v, q, m = part.x, part.v, part.q, part.m
  F = q*E # Lorentz force, F=qE
  a = F/m # acceleration, F=ma

  v[1:np,:] .= v[1:np,:] .+ Δt * a[1:np,:]
  x[1:np,:] .= x[1:np,:] .+ Δt * v[1:np,:]

  Diagnostics.register_diagnostic("pv"*part.name, ParticleVectorData(x[1:np,1], x[1:np,2], v[1:np,1], v[1:np,2], zeros(np), id[1:np], 1))
  Diagnostics.register_diagnostic("pE"*part.name, ParticleVectorData(x[1:np,1], x[1:np,2], E[1:np,1], E[1:np,2], zeros(np), id[1:np], 1))
end
end
