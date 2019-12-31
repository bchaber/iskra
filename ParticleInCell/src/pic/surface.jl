using DataStructures

Base.findfirst(x::Array{Int64,1}, A::Array{Int64,2}) =
  for i=1:size(A,2)
    if x == A[:,i]
    return i  
    end
  end
Base.in(x::Array{Int64,1}, A::Array{Int64,2}) =
  findfirst(x, A) ≠ nothing

abstract type Surface end
mutable struct ReflectiveSurface <: Surface end
mutable struct AbsorbingSurface <: Surface end

create_absorbing_surface()  = AbsorbingSurface()
create_reflective_surface() = ReflectiveSurface()
function absorbs(s::Surface)           false end
function absorbs(s::ReflectiveSurface) false end
function absorbs(s::AbsorbingSurface)  true  end
function hit!(s::Surface, part::KineticSpecies, p::Int64,
              Δt::Float64, n̂::Array{Float64,1})
end
function hit!(s::ReflectiveSurface, part::KineticSpecies, p::Int64,
              Δt::Float64, n̂::Array{Float64,1})
  part.x[p,:] .-= part.v[p,:]*Δt
  for i=[1,2,3]
    if n̂[i] ≠ 0
      part.v[p,i] *= -1
    end
  end
  part.x[p,:] .+= part.v[p,:]*Δt
end
function hit!(s::AbsorbingSurface, part::KineticSpecies, p::Int64,
              Δt::Float64, n̂::Array{Float64,1})
end

struct SurfaceTracker
  cells :: Array{Int64,2}
  tracked :: Array{Tuple{Int64,Int64,Int64,Float64,Float64,Float64},1}
  surfaces :: Array{Surface,1}
  boundaries :: Array{Int64,1}
  Δh :: Float64
  Δt :: Float64
end

function create_surface_tracker(bcs::Array{Int8, 3}, surfaces::Array{Surface,1}, Δh, Δt)
  nx, ny = size(bcs)
  cells = Array{Int64,1}[]
  tracked = Tuple{Int64,Int64,Int64,Float64,Float64,Float64}[]
  boundaries  = Int64[]
  
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

      if a != 0
        if a == b
          push!(cells, [i,j, i,  j-1])
          push!(boundaries,  a)
        else
          push!(cells, [i,j, i-1,j-1])
          push!(boundaries,  a)
        end
      end

      if b != 0
        if b == c
          push!(cells, [i,j, i+1,j])
          push!(boundaries,  b)
        else
          push!(cells, [i,j, i+1,j-1])
          push!(boundaries,  b)
        end
      end

      if c != 0
        if c == d
          push!(cells, [i,j, i,  j+1])
          push!(boundaries,  c)
        else
          push!(cells, [i,j, i+1,j+1])
          push!(boundaries,  c)
        end
      end

      if d != 0
        if d == a
          push!(cells, [i,j, i-1,j])
          push!(boundaries,  d)
        else
          push!(cells, [i,j, i-1,j+1])
          push!(boundaries,  d)
        end
      end
    end
  end
  SurfaceTracker(hcat(cells...), tracked, surfaces, boundaries, Δh, Δt)
end

function track!(::Nothing, part::KineticSpecies) end
function track!(st::SurfaceTracker, part::KineticSpecies)
  Δh = st.Δh
  Δt = st.Δt
  np = part.np
  px = view(part.x, 1:np, :)
  pv = view(part.v, 1:np, :)

  deleteat!(st.tracked, 1:length(st.tracked))
  for p=1:part.np
    i, j, hx, hy = particle_cell(px, p, Δh)
    if [i,j] ∈ st.cells[1:2,:] ||
       [i,j] ∈ st.cells[3:4,:]
      push!(st.tracked, (p,i,j,hx,hy,Δt))
    end
  end
end

function check_particle(i, j, hx, hy, Δt, Δh, vx, vy)
  dx = (vx > 0) ? Δh*(1-hx) : Δh*hx
  dy = (vy > 0) ? Δh*(1-hy) : Δh*hy

  Δtx, Δty = dx/abs(vx), dy/abs(vy)
  if Δt < Δtx && Δt < Δty # particle stayed in the same cell
    return i, j, hx, hy, Δt
  end

  if Δtx < Δty
    return (vx > 0) ?
      (i+1, j, 0., hy + vy*Δtx/Δh, Δt-Δtx) :
      (i-1, j, 1., hy + vy*Δtx/Δh, Δt-Δtx)
  else
    return (vy > 0) ?
      (i, j+1, hx + vx*Δty/Δh, 0., Δt-Δty) :
      (i, j-1, hx + vx*Δty/Δh, 1., Δt-Δty)
  end
end
    
function check!(::Nothing, part::KineticSpecies, Δt) end
function check!(st::SurfaceTracker, part::KineticSpecies, Δt)
  Δh = st.Δh
  px = view(part.x, 1:part.np, :)
  pv = view(part.v, 1:part.np, :)
  absorbed = SortedSet{Int64}(Base.Order.Reverse)

  while length(st.tracked) > 0
    p, i , j , hx , hy , Δt  = popfirst!(st.tracked)
    vx, vy = pv[p,1], pv[p,2]

    ΔTx, ΔTy = Δh/abs(vx), Δh/abs(vy)
    if ΔTx < Δt || ΔTx < Δt
      println("ERROR: Particle is too fast!")
    end
    
    i′, j′, hx′, hy′, Δt′ = check_particle(i, j, hx, hy, Δt, Δh, vx, vy)
    if i == i′ && j == j′
      continue
    end
    n̂ = [i′-i, j′-j, 0.]
    r = findfirst([i, j, i′,j′], st.cells)
    r = (r ≠ nothing) ? r :
        findfirst([i′,j′,i, j],  st.cells)

    if r ≠ nothing
      boundary = st.boundaries[r] + 1
      surface = st.surfaces[boundary]
      hit!(surface, part, p, Δt′, n̂)
      if absorbs(surface)
        push!(absorbed, p)
      end
    else
      push!(st.tracked, (p, i′, j′, hx′, hy′, Δt′))
    end
  end

  for p in absorbed
    remove!(part, p)
  end
end