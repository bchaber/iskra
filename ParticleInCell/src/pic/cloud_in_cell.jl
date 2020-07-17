function particle_to_grid(part::KineticSpecies{1, V}, grid::CartesianGrid{1}, pu) where V
  nx, = size(grid)
  Δx, = grid.Δh
  np = part.np
  px = view(part.x, 1:np, :)
  u  = zeros(nx)   # distribution of pu on grid

  for p=1:np
    vu = view(pu, p)
    (i,), (hx,) = particle_cell(px, p, grid.Δh)
    # interpolate charge to nodes
    u[i  ,1] += (1-hx)*vu/Δx
    u[i+1,1] +=    hx *vu/Δx
  end

  return u
end

function particle_to_grid(part::KineticSpecies{2, V}, grid::CartesianGrid{2}, pu) where V
  nx, ny = size(grid)
  Δx, Δy = grid.Δh
  np = part.np
  px = view(part.x, 1:np, :)
  u  = zeros(nx, ny)   # distribution of pu on grid

  for p=1:np
    vu = view(pu, p)
    (i, j), (hx, hy) = particle_cell(px, p, grid.Δh)
    # interpolate charge to nodes
    u[i  ,j]   += (1-hx)*(1-hy)* vu/(Δx*Δy)
    u[i+1,j]   +=    hx *(1-hy)* vu/(Δx*Δy)
    u[i  ,j+1] += (1-hx)*   hy * vu/(Δx*Δy)
    u[i+1,j+1] +=    hx *   hy * vu/(Δx*Δy)
  end

  return u
end

function grid_to_particle(grid::CartesianGrid{1}, part::KineticSpecies{1, V}, u) where V
  np = part.np
  Δh = grid.Δh
  px = view(part.x, 1:np, :)
  pu = zeros(np, V)

  for p=1:np
    (i,), (hx,) = particle_cell(px, p, Δh)
    # gather electric field
    pu[p,:] = (1.0-hx)*view(u, i,   :) +
              (1.0-hx)*view(u, i+1, :)
  end

  return pu
end

function grid_to_particle(grid::CartesianGrid{2}, part::KineticSpecies{2, V}, u) where V
  np = part.np
  Δh = grid.Δh
  px = view(part.x, 1:np, :)
  pu = zeros(np, V)

  for p=1:np
    (i, j), (hx, hy) = particle_cell(px, p, Δh)
    # gather electric field
    pu[p,:] = (1.0-hx)*(1.0-hy)* view(u, i , j,  :) + # contribution from (i,j)
                   hx *(1.0-hy)* view(u, i+1,j,  :) + # (i+1,j)
              (1.0-hx)*     hy * view(u, i,  j+1,:) + # (i,j+1)
                   hx *     hy * view(u, i+1,j+1,:)   # (i+1,j+1)
  end

  return pu
end
