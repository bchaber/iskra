

function wrap!(part::KineticSpecies, grid::CartesianGrid{2})
  nx, ny = grid.n .- 1
  Δx, Δy = grid.Δh
  Lx, Ly = [nx, ny] .* [Δx, Δy]
  ox, oy = grid.origin
  bb = [ox oy] .+ [0. Lx; 0. Ly]
  px = view(part.x, 1:part.np, :)
  for p=1:part.np
    if px[p,1] < bb[1,1] px[p,1] += Lx end
    if px[p,1] > bb[1,2] px[p,1] -= Lx end
    if px[p,2] < bb[2,1] px[p,2] += Ly end
    if px[p,2] > bb[2,2] px[p,2] -= Ly end
  end
end