function particle_to_grid!(u, part::KineticSpecies{1, V}, grid::CartesianGrid{1}, pu) where V
  for p=1:part.np
    i, _, hx, _ = particle_cell(part, p, grid.Δh)
    u[i] += pu[p]
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
