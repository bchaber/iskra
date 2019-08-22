function particle_to_grid(part, grid, pu)
  nx, ny = size(grid)
  Δh = grid.Δh
  np = part.np
  px = view(part.x, 1:np, :)
  u  = zeros(nx, ny)   # distribution of pu on grid

  for p=1:np
    i, j, hx, hy = particle_cell(px, p, Δh)
    # interpolate charge to nodes
    u[i  ,j]   += (1-hx)*(1-hy)* pu(p)/Δh^2
    u[i+1,j]   +=    hx *(1-hy)* pu(p)/Δh^2
    u[i  ,j+1] += (1-hx)*   hy * pu(p)/Δh^2
    u[i+1,j+1] +=    hx *   hy * pu(p)/Δh^2
  end

  u[ 1,:] .*= 2
  u[nx,:] .*= 2
  u[:, 1] .*= 2
  u[:,ny] .*= 2

  return u
end

function grid_to_particle(grid, part, u)
  np = part.np
  Δh = grid.Δh
  px = view(part.x, 1:np, :)
  pu = zeros(np, 3)

  for p=1:np
    i, j, hx, hy = particle_cell(px, p, Δh)
    # gather electric field
    pu[p,:] = (1-hx)*(1-hy)* u(i , j  ) + # contribution from (i,j)
                 hx *(1-hy)* u(i+1,j  ) + # (i+1,j)
              (1-hx)*   hy * u(i,  j+1) + # (i,j+1)
                 hx *   hy * u(i+1,j+1)   # (i+1,j+1)
  end

  return pu
end
