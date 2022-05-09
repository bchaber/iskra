#        .
#       / \
#      /   \
#     /     \
# -+-o-+-o-+-o-+-

function particle_to_grid!(u, part::KineticSpecies{1, V}, grid::CartesianGrid{1}, pu) where V
  n = grid.n[1] - 1
  for p=1:part.np
    i, _, hx, _ = particle_cell(part, p, grid.Δh)
    b = abs(0.5 - hx)
    a = 1.0 - b

    if i == 1 || i == n
      u[i] += hx > 0.5 ? b * pu[p] : a * pu[p]
    end
    
    if hx > 0.5
      u[i]   += b * pu[p]
      u[i+1] += a * pu[p]
    else
      u[i-1] += b * pu[p]
      u[i]   += a * pu[p]
    end
  end

  return u
end

function grid_to_particle!(pu::Vector{SVector{V, Float64}}, grid::CartesianGrid{1}, part::KineticSpecies{1, V}, u) where V
  n = grid.n[1] - 1
  for p=1:part.np
    i, _, hx, _ = particle_cell(part, p, grid.Δh)
    b = abs(0.5 - hx)
    a = 1.0 - b

    if i == 1 || i == n
      pu[p] += hx > 0.5 ? b * u[i] : a * u[i]
    end
    
    if hx > 0.5
      pu[p] += b * u[i]
      pu[p] += a * u[i+1]
    else
      pu[p] += b * u[i-1]
      pu[p] += a * u[i]
  end

  return pu
end