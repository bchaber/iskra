

function wrap!(part::KineticSpecies, grid::UniformGrid{XY2D})
  nx, ny, nz = grid.nx-1, grid.ny-1, grid.nz-1
  Δx, Δy, Δz = grid.Δh
  Lx, Ly, Lz = [nx, ny, nz] .* [Δx, Δy, Δz]
  origin = grid.origin
  bb = [origin origin] .+ [0. Lx; 0. Ly; 0 Lz]
  px = view(part.x, 1:part.np, :)
  for p=1:part.np
    if px[p,1] < bb[1,1] px[p,1] += Lx end
    if px[p,1] > bb[1,2] px[p,1] -= Lx end
    if px[p,2] < bb[2,1] px[p,2] += Ly end
    if px[p,2] > bb[2,2] px[p,2] -= Ly end
  end
end