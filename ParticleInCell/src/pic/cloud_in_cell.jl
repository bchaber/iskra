function particle_to_grid!(u, part::KineticSpecies{1, V}, grid::CartesianGrid{1}, pu) where V
  n = grid.n[1] - 1
  for p=1:part.np
    i, _, hx, _ = particle_cell(part, p, grid.Δh)
    if hx == 0.5
      u[i] += pu[p]
    elseif hx < 0.5
      u[i] += (0.5 + hx) * pu[p]
      if i > 1
        u[i-1] += (0.5 - hx) * pu[p]
      else
        u[i]   += (0.5 - hx) * pu[p]
      end
    elseif hx > 0.5
      u[i] += (1.5 - hx) * pu[p]
      if i < n
        u[i+1] += (hx - 0.5) * pu[p]
      else
        u[i]   += (hx - 0.5) * pu[p]
      end
    end
  end

  return u
end

function grid_to_particle!(pu::Vector{SVector{V, Float64}}, grid::CartesianGrid{1}, part::KineticSpecies{1, V}, u) where V
  for p=1:part.np
    i, _, hx, _ = particle_cell(part, p, grid.Δh)
    pu[p] = (1.0-hx) * u[i] + (hx) * u[i+1]
  end

  return pu
end
