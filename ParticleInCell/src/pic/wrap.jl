function discard!(part::KineticSpecies{D,V}, grid::UniformGrid{CS, D}, dims) where {CS, D, V}
  before = part.np
  for i=dims
    nx = grid.n[i] - 1
    Δx = grid.Δh[i]
    Lx = nx*Δx
    ox = grid.origin[i]
    px = part.x
    for p=reverse(1:part.np)
      α = fld(px[p][i] - ox, Lx)
      if α ≠ 0
        remove!(part, p)
      end
    end
  end
  after = part.np
  return before - after
end

function wrap!(part::KineticSpecies{D,V}, grid::UniformGrid{CS, D}, dx, dims) where {CS, D, V}
  for i=dims
    nx = grid.n[i] - 1
    Δx = grid.Δh[i]
    Lx = nx*Δx
    ox = grid.origin[i]
    px = part.x
    for p=1:part.np
      α = fld(px[p][i] - ox, Lx)
      #if p == 1190
      #  println(part, " wrap: ", α, " x: ", px[p], " ", α * Lx * dx)
      #end
      if α ≠ 0
        px[p] -= α * Lx * dx
      end
    end
  end
end