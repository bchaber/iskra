function discard!(part::KineticSpecies{D,V}, grid::UniformGrid{CS, D}; dims=1:D) where {CS, D, V}
  before = part.np
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
  after = part.np
  return before - after
end

function wrap!(part::KineticSpecies{D,V}, grid::UniformGrid{CS, D}; dims=1:D) where {CS, D, V}
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