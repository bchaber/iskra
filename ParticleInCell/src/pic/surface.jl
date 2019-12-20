Base.findfirst(x::Array{Int64,1}, A::Array{Int64,2}) =
  for i=1:size(A,2)
    if x == A[:,i]
    return i  
    end
  end
Base.in(x::Array{Int64,1}, A::Array{Int64,2}) =
  findfirst(x, A) ≠ nothing

abstract type Surface end
mutable struct MetalSurface <: Surface
  hits :: Int64
end
mutable struct AbsorbingSurface <: Surface
  hits :: Int64
end
struct SurfaceTracker
  cells :: Array{Int64,2}
  boundaries :: Array{Int64,1}
  tracked :: Array{Tuple{Int64,Int64,Int64},1}
  Δh :: Float64
end

function create_surface_tracker(bcs::Array{Int8, 3}, Δh)
  nx, ny = size(bcs)
  cells = Array{Int64,1}[]
  boundaries  = Int64[]
  tracked = Tuple{Int64,Int64,Int64}[]
  
  for i=1:nx-1
    for j=1:ny-1
      a = bcs[i,  j  ,1]
      b = bcs[i+1,j  ,1]
      c = bcs[i+1,j+1,1]
      d = bcs[i  ,j+1,1]
      if j == 1
        push!(cells, [i,j-1, i,j])
        push!(boundaries,  (a == b) ? a : 0)
      end

      if i == nx-1
        push!(cells, [i,j, i+1,j])
        push!(boundaries,  (b == c) ? b : 0)
      end

      if j == ny-1
        push!(cells, [i,j, i,j+1])
        push!(boundaries,  (c == d) ? c : 0)
      end

      if i == 1
        push!(cells, [i-1,j, i,j])
        push!(boundaries,  (d == a) ? d : 0)
      end
    end
  end
  println("Created tracker for boundaries")
  for k=1:length(boundaries)
    println(cells[k], " ", boundaries[k])
  end
  SurfaceTracker(hcat(cells...), boundaries, tracked, Δh)
end

function track!(::Nothing, part::KineticSpecies) end
function track!(st::SurfaceTracker, part::KineticSpecies)
  Δh = st.Δh
  np = part.np
  px = view(part.x, 1:np, :)
  pv = view(part.v, 1:np, :)

  deleteat!(st.tracked, 1:length(st.tracked))
  for p=1:part.np
    i, j, hx, hy = particle_cell(px, p, Δh)
    if [i,j] ∈ st.cells[1:2,:] ||
       [i,j] ∈ st.cells[3:4,:]
      push!(st.tracked, (p,i,j))
      println("Started tracking (", p ,"|", i, ",", j, ") ", hx, " ", hy)
    end
  end
end

function check!(::Nothing, part::KineticSpecies, Δt) end
function check!(st::SurfaceTracker, part::KineticSpecies, Δt)
  Δh = st.Δh
  px = view(part.x, 1:part.np, :)
  pv = view(part.v, 1:part.np, :)
  
  for (p,i,j) in st.tracked
    ~, ~, hx, hy = particle_cell(px, p, Δh)
    vx, vy = pv[p,1], pv[p,2]
    dx = (vx > 0) ? Δh*(1-hx) : Δh*hx
    dy = (vy > 0) ? Δh*(1-hy) : Δh*hy

    ΔTx, ΔTy = Δh/abs(vx), Δh/abs(vy)
    if ΔTx < Δt || ΔTx < Δt
      println("ERROR: Particle ", p, " is too fast!")
    end

    Δtx, Δty = dx/abs(vx), dy/abs(vy)
    println("Δtx ", Δtx/Δt, " Δty ", Δty/Δt)
    if Δt < Δtx && Δt < Δty # particle stayed in the same cell
      continue
    end

    cells = (Δtx < Δty) ?
      (vx < 0 ? [i-1,j,i,j] : [i,j,i+1,j]) :
      (vy < 0 ? [i,j-1,i,j] : [i,j,i,j+1])
    println("Looking for ", cells)
    r = findfirst(cells, st.cells)
    if r ≠ nothing
      println("Surface ", st.boundaries[r], " hit!")
      remove!(part, p)
    end
  end
end