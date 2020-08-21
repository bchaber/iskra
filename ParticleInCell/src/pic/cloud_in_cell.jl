function particle_to_grid(part::KineticSpecies{2, V}, grid::CartesianGrid{2}, pu) where V
  nx, ny = size(grid)
  Δx, Δy = grid.Δh
  np = part.np
  px = view(part.x, 1:np, :)
  u  = zeros(nx, ny)   # distribution of pu on grid

  for p=1:np
    i, j, hx, hy = particle_cell(px, p, grid.Δh)
    # interpolate charge to nodes
    u[i  ,j]   += (1.0-hx) * (1.0-hy) * pu(p)
    u[i+1,j]   +=     (hx) * (1.0-hy) * pu(p)
    u[i  ,j+1] += (1.0-hx) *     (hy) * pu(p)
    u[i+1,j+1] +=     (hx) *     (hy) * pu(p)
  end

  return u
end

function grid_to_particle(grid::CartesianGrid{2}, part::KineticSpecies{2, V}, u) where V
  np = part.np
  Δh = grid.Δh
  px = view(part.x, 1:np, :)
  pu = zeros(np, V)

  for p=1:np
    i, j, hx, hy = particle_cell(px, p, Δh)
    # gather electric field
    pu[p,:] = (1.0-hx) * (1.0-hy) * u(i , j  ) + # contribution from (i,j)
                  (hx) * (1.0-hy) * u(i+1,j  ) + # (i+1,j)
              (1.0-hx) *     (hy) * u(i,  j+1) + # (i,j+1)
                  (hx) *     (hy) * u(i+1,j+1)   # (i+1,j+1)
  end

  return pu
end

function particle_to_grid(part::KineticSpecies{2, V}, grid::AxialGrid{2}, pu) where V
  nr, nz = size(grid)
  Δr, Δz = grid.Δh
  np = part.np
  px = view(part.x, 1:np, :)
  u  = zeros(nr, nz)   # distribution of pu on grid

  for p=1:np
    i, j, hr, hz = particle_cell(px, p, grid.Δh)
    # interpolate charge to nodes
    u[i  ,j]   += (1.0-hr) * (1.0-hz) * pu(p)
    u[i+1,j]   +=     (hr) * (1.0-hz) * pu(p)
    u[i  ,j+1] += (1.0-hr) *     (hz) * pu(p)
    u[i+1,j+1] +=     (hr) *     (hz) * pu(p)
  end

  return u
end

function grid_to_particle(grid::AxialGrid{2}, part::KineticSpecies{2, V}, u) where V
  np = part.np
  px = view(part.x, 1:np, :)
  pu = zeros(np, V)

  Δr, Δz = grid.Δh
  for p=1:np
    i, j, hr, hz = particle_cell(px, p, grid.Δh)
    # gather electric field
    pu[p,:] = (1.0-hr) * (1.0-hz) * u(i , j  ) + # contribution from (i,j)
                  (hr) * (1.0-hz) * u(i+1,j  ) + # (i+1,j)
              (1.0-hr) *     (hz) * u(i,  j+1) + # (i,j+1)
                  (hr) *     (hz) * u(i+1,j+1)   # (i+1,j+1)
  end

  return pu
end
