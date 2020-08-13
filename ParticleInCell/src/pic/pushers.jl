export create_boris_pusher, create_axial_boris_pusher

struct BorisPusher{T} end
create_boris_pusher() = BorisPusher{:xy}()
create_axial_boris_pusher() = BorisPusher{:rz}()

function push_particles!(::BorisPusher{:xy},
	part::KineticSpecies{D, V}, E, Δt) where {D, V}
  push_in_cartesian!(part, E, Δt)
end

function push_particles!(::BorisPusher{:rz},
  part::KineticSpecies{D, V}, E, Δt) where {D, V}
  push_in_cartesian!(part, E, Δt)
  transform_from_cartesian_to_cylindrical!(part, Δt)
end

function push_in_cartesian!(part::KineticSpecies{D, V}, E, Δt) where {D, V}
  np, id = part.np, part.id
  x, v, q, m = part.x, part.v, part.q, part.m
  F = q*E # Lorentz force, F=qE
  a = F/m # acceleration,  F=ma

  v[1:np,1:V] .= @views v[1:np,1:V] .+ Δt * a[1:np,1:V]
  x[1:np,1:D] .= @views x[1:np,1:D] .+ Δt * v[1:np,1:D]
end

function transform_from_cartesian_to_cylindrical!(part::KineticSpecies{2, 3}, Δt)
  x, v, np = part.x, part.v, part.np
  y = Δt * v[1:np,3] # y ← Δt⋅vz
  r = sqrt.(x[1:np,1].^2 .+ y[1:np].^2) # r ← √(x² + y²)

  sinθ = y[1:np,1]./r
  cosθ = sqrt.(1.0 .- sinθ.^2)
  vr = cosθ .* v[1:np,1] .+ sinθ .* v[1:np,3]
  vy =-sinθ .* v[1:np,1] .+ cosθ .* v[1:np,3]
  x[1:np,1] .= r
  v[1:np,1] .= vr
  v[1:np,3] .= vy
end