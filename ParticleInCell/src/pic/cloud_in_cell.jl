#        .
#       / \
#      /   \
#     /     \
# -+-o-+-o-+-o-+-

function particle_to_grid!(u, part::KineticSpecies{1, V}, grid::CartesianGrid{1}, pu) where V
  nx = grid.n[1] - 1
  for p=1:part.np
    i, _, hx, _ = particle_cell(part, p, grid.Δh)
    b = abs(0.5 - hx)
    a = 1.0 - b

    if hx > 0.5
      u[i]   += a * pu[p]
      u[i < nx ? i+1 : 1] += b * pu[p]
    else
      u[i > 1 ? i-1 : nx] += b * pu[p]
      u[i]   += a * pu[p]
    end
  end

  return u
end

function grid_to_particle!(pu::Vector{SVector{V, Float64}}, grid::CartesianGrid{1}, part::KineticSpecies{1, V}, u) where V
  nx = grid.n[1] - 1
  for p=1:part.np
    i, _, hx, _ = particle_cell(part, p, grid.Δh)
    b = abs(0.5 - hx)
    a = 1.0 - b

    if hx > 0.5
      pu[p] = b * u[i < nx ? i+1 : 1] + a * u[i]
    else
      pu[p] = b * u[i > 1 ? i-1 : nx] + a * u[i]
    end
  end

  return pu
end