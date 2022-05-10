struct BorisPusherData{T}
  E :: T # electric field
  B :: T # magnetic field
  V :: T # electric field integral
  I :: T # magnetic field integral
end

struct BorisPusher{CS}
  data :: Dict{Symbol, BorisPusherData{Vector{SVector{3, Float64}}}}
end

function create_boris_pusher(species)
  data = Dict{Symbol,BorisPusherData}()
  for part in species
    if is_fluid(part) continue end
    data[part.name] = BorisPusherData(
      similar(part.v), similar(part.v),
      similar(part.v), similar(part.v),
    )
  end
  BorisPusher{:xy}(data)
end

function push_particles!(pusher::BorisPusher{:xy},
	part::KineticSpecies{D, V}, Δt) where {D, V}
  data = pusher.data[part.name]
  push_in_cartesian!(part, data, Δt)
end

function push_in_cartesian!(part::KineticSpecies{1,3}, data::BorisPusherData, Δt)
  E, B = data.E, data.B
  V, I = data.V, data.I
  x, v = part.x, part.v
  for i in 1:part.np
    V[i] = 0.5Δt * part.q/part.m * E[i]
    I[i] = 0.5Δt * part.q/part.m * B[i]

    b  = B[i]
    t  = I[i]
    v⁻ = V[i] + v[i]
    v′ = v⁻ + v⁻ × b

    s  = 2.0 / (1.0 + (t ⋅ t)) * t
    v⁺ = v⁻ + v′ × s
    #if any(isnan.(V[i] + v⁺))
    #  println("E: ", E[i], " B: ", B[i])
    #  println("V: ", V[i], " V: ", I[i])
    #  println("b: ", b, " t: ", t)
    #  println("v⁻: ", v⁻, " v′: ", v′)
    #end

    #if i == 1190
    #  println(part, " vel: ", v[i][1], " => ", V[i][1] + v⁺[1])
    #  println(part, " pos: ", x[i][1], " => ", x[i][1] + Δt * (V[i][1] + v⁺[1]))
    #end
    v[i] = V[i] + v⁺
    x[i] += Δt * @SVector[v[i][1]]
  end
  return nothing
end
