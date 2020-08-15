export create_boris_pusher, create_axial_boris_pusher
export initial_half_step

struct BorisPusher{T} end
create_boris_pusher() = BorisPusher{:xy}()
create_axial_boris_pusher() = BorisPusher{:rz}()

function push_particles!(::BorisPusher{:xy},
	part::KineticSpecies{D, V}, E, B, Δt) where {D, V}
  push_in_cartesian!(part, E, B, Δt)
end

function push_particles!(::BorisPusher{:rz},
  part::KineticSpecies{D, V}, E, B, Δt) where {D, V}
  push_in_cartesian!(part, E, B, Δt)
  transform_from_cartesian_to_cylindrical!(part, Δt)
end

function LinearAlgebra.dot(u::Matrix{Float64}, v::Matrix{Float64})
  n, m = size(u)
  w = zeros(n)
  for i=1:n
    w[i,:] .= @views u[i,:] ⋅ v[i,:]
  end
  return w
end

function LinearAlgebra.cross(u::Matrix{Float64}, v::Matrix{Float64})
  n, m = size(u)
  w = zeros(n, m)
  for i=1:n
    w[i,:] .= @views u[i,:] × v[i,:]
  end
  return w
end

function push_in_cartesian!(part::KineticSpecies{D, V}, E, B, Δt) where {D, V}
  np, id = part.np, part.id
  qm = part.q/part.m
  x, v = view(part.x, 1:np, :), view(part.v, 1:np, :)
  v⁻ = 0.5Δt * qm * E .+ v
  t  = 0.5Δt * qm * B
  t² = t ⋅ t
  v′ = v⁻ .+ v⁻ × B
  s  = 2.0 ./ (1.0 .+ t²) .* t
  v⁺ = v⁻ .+ v′ × s

  v[:,1:V] .= Δt * E * qm * 0.5 .+ v⁺
  x[:,1:D] .= Δt * v[:,1:D] .+ x[:,1:D]
end

function transform_from_cartesian_to_cylindrical!(part::KineticSpecies{2, 3}, Δt)
  x, v, np = part.x, part.v, part.np
  y = Δt * v[1:np,3] # y ← Δt⋅vz
  r = sqrt.(x[1:np,1].^2 .+ y[1:np].^2) # r ← √(x² + y²)

  sinθ = y[1:np,1]./r
  sinθ[r .≈ 0.0] .= 0.0
  cosθ = sqrt.(1.0 .- sinθ.^2)

  vr = cosθ .* v[1:np,1] .+ sinθ .* v[1:np,3]
  vy =-sinθ .* v[1:np,1] .+ cosθ .* v[1:np,3]
  x[1:np,1] .= r
  v[1:np,1] .= vr
  v[1:np,3] .= vy
end

function initial_half_step(v0, qm, E, Δt)
  v = v0 - qm * E * 0.5Δt
end