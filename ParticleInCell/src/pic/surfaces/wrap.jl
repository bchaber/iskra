function discard!(part::KineticSpecies{D,V}, grid::CartesianGrid{D}; dims=1:D) where {D, V}
  for i=dims
    nx = grid.n[i] - 1
    Δx = grid.Δh[i]
    Lx = nx*Δx
    ox = grid.origin[i]
    px = view(part.x, 1:part.np, i)
    for p=reverse(1:part.np)
      α = fld(px[p] - ox, Lx)
      if α ≠ 0
        remove!(part, p)
      end
    end
  end
end

function wrap!(part::KineticSpecies{D,V}, grid::CartesianGrid{D}; dims=1:D) where {D, V}
  for i=dims
    nx = grid.n[i] - 1
    Δx = grid.Δh[i]
    Lx = nx*Δx
    ox = grid.origin[i]
    px = view(part.x, 1:part.np, i)
    for p=1:part.np
      α = fld(px[p] - ox, Lx)
      if α ≠ 0
        px[p] -= α*Lx
      end
    end
  end
end