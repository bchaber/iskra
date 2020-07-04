using RegularGrid

mutable struct DensitySource{D}
  δ :: Array{Float64, D}
  grid :: CartesianGrid{D}
end

mutable struct MaxwellianSource{D, V}
   rate :: Float64
   x :: AbstractArray{Float64, 2}
   v :: AbstractArray{Float64, 2}
   species
end
MaxwellianSource(rate, x, v) =
MaxwellianSource{2,3}(float(rate), float.(x), v, nothing)
  
function sample!(config :: MaxwellianSource{2,3}, species :: KineticSpecies{2,3}, Δt)
  np = species.np
  cx, cv = config.x, config.v
  px, pv = @views species.x[1+np:end,:], species.v[1+np:end,:]

  n = minimum([size(px, 1), floor(Integer, config.rate*Δt)])
  px[1:n,1]=rand(n,1)*cx[1,2]  .+ cx[1,1]
  px[1:n,2]=rand(n,1)*cx[2,2]  .+ cx[2,1]
  pv[1:n,1]=randn(n,1)*cv[1,2] .+ cv[1,1]
  pv[1:n,2]=randn(n,1)*cv[2,2] .+ cv[2,1]

  species.np += n
end

function sample!(config :: DensitySource{D}, species :: FluidSpecies{D}, Δt) where D
  species.n .+= config.δ
end

function create!(species :: KineticSpecies, grid, n)
  np = species.np
  px, pv, pw = @views species.x[1+np:end,:], species.v[1+np:end,:], species.wg[1+np:end]
  nx, ny = size(n)
  Δx, Δy, ~ = grid.Δh
  N = ceil.(max.(n/species.w0, 0))
  N = ceil(Int64, sum(N))
  pv[1:N,1]=rand(N,1).+rand(N,1).+rand(N,1).-1.5
  pv[1:N,2]=rand(N,1).+rand(N,1).+rand(N,1).-1.5
  px[1:N,1]=rand(N,1)*Δx.-Δx/2 # relative x position in cell
  px[1:N,2]=rand(N,1)*Δy.-Δy/2 # relative y position in cell
  s = 1
  for i=2:nx-1
    for j=2:ny-1
      if n[i,j] > species.w0
        dn = floor(Int64, n[i,j]/species.w0)
        px[s:s+dn-1,1] .+= grid.x[i,j]
        px[s:s+dn-1,2] .+= grid.y[i,j]
        n[i,j] -= dn*species.w0
        s += dn
      end
      if n[i,j] > 1e-6species.w0
        px[s,1] += grid.x[i,j]
        px[s,2] += grid.y[i,j]
        pw[s] = n[i,j]
        n[i,j] = 0
        s += 1
      end
    end
  end
  species.np += (s-1)
end

function destroy!(species :: KineticSpecies, grid, n)
  np, p = species.np, 1
  px, pw = @views species.x[1:np,:], species.wg[1:np]
  while p <= np
    i, j, _, _ = particle_cell(px, p, grid.Δh)
    if n[i,j] >= 0
      p = p + 1
      continue
    end

    if pw[p] < -n[i,j]
      n[i,j] += pw[p]
      remove!(species, p)
      np -= 1
      continue
    end

    if pw[p] > -n[i,j]
      pw[p] -= n[i,j]
      n[i,j] = 0
    end
    p = p + 1
  end
end