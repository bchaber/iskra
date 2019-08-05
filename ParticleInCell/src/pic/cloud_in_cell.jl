function particle_to_grid(part, grid)
  nx, ny = size(grid)
  np = part.np
  Δh = grid.Δh
  px = part.x[1:np,:]
  ρ  = zeros(nx, ny)   # charge distribution

  for p=1:np
    i, j, hx, hy = particle_cell(px, p, Δh)
    np2c, q = part.np2c, part.q
    # interpolate charge to nodes
    ρ[i  ,j]   += (1-hx)*(1-hy)* np2c * q / Δh^2
    ρ[i+1,j]   +=    hx *(1-hy)* np2c * q / Δh^2
    ρ[i  ,j+1] += (1-hx)*   hy * np2c * q / Δh^2
    ρ[i+1,j+1] +=    hx *   hy * np2c * q / Δh^2
  end

  ρ[1,:]  = 2ρ[1,:]
  ρ[nx,:] = 2ρ[nx,:]
  ρ[:,1]  = 2ρ[:,1]
  ρ[:,ny] = 2ρ[:,ny]

  return ρ
end

function grid_to_particle(grid, part, Ex, Ey)
  np = part.np
  Δh = grid.Δh
  px = part.x[1:np,:]
  pE = zeros(np, 2)

  for p=1:np
    i, j, hx, hy = particle_cell(px, p, Δh)
    # gather electric field
    pE[p,:] = (1-hx)*(1-hy)*[Ex[i,  j  ] Ey[i,  j  ]] + # contribution from (i,j)
                 hx *(1-hy)*[Ex[i+1,j  ] Ey[i+1,j  ]] + # (i+1,j)
              (1-hx)*   hy *[Ex[i,  j+1] Ey[i,  j+1]] + # (i,j+1)
                 hx *   hy *[Ex[i+1,j+1] Ey[i+1,j+1]]   # (i+1,j+1)
  end
  return pE
end
