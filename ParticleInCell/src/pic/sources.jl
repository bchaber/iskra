using RegularGrid

mutable struct DensitySource
  δ :: Array{Float64,2}
  grid :: UniformGrid
end

mutable struct MaxwellianSource
   species
   rate :: Float64
   x :: AbstractArray{Float64,2}
   v :: AbstractArray{Float64,2}
   n :: Float64
end

function create_maxwellian_source(species, rate, x, v)
  MaxwellianSource(species, rate, x, v, 0.0)
end
  
function sample!(config :: MaxwellianSource, species :: KineticSpecies, Δt)
  np = species.np
  cx, cv = config.x, config.v
  px, pv = @views species.x[1+np:end,:], species.v[1+np:end,:]

  config.n += config.rate*Δt
  n = minimum([size(px, 1), floor(Integer, config.n)])
  px[1:n,1]=randn(n,1)*cx[1,2]/10 .+ cx[1,2]/2 # x position
  px[1:n,2]=randn(n,1)*cx[2,2]/10 .+ cx[2,2]/2 # y position
  pv[1:n,1]=1.0*(-1.5.+rand(n,1).+rand(n,1).+rand(n,1))*cv[1,1] .+ cv[1,2]
  pv[1:n,2]=0.5*(-1.5.+rand(n,1).+rand(n,1).+rand(n,1))*cv[2,1] .+ cv[2,2]

  config.n   -= n
  species.np += n
end

function sample!(config :: DensitySource, species :: FluidSpecies, Δt)
  species.n .+= config.δ
end

function sample!(config :: DensitySource, species :: KineticSpecies, Δt)
  np = species.np
  px, pv = @views species.x[1+np:end,:], species.v[1+np:end,:]
  δ = round.(Integer, config.δ * config.grid.Δh^2/species.np2c)
  n = minimum([size(px, 1), sum(max.(δ, 0))])
  px[1:n,1]=rand(n,1)*config.grid.Δh # relative x position in cell
  px[1:n,2]=rand(n,1)*config.grid.Δh # relative y position in cell
  pv[1:n,1]=1.0*(-1.5.+rand(n,1).+rand(n,1).+rand(n,1))
  pv[1:n,2]=0.5*(-1.5.+rand(n,1).+rand(n,1).+rand(n,1))
  s = 1
  nx, ny = size(config.δ)
  for i=1:nx
    for j=1:ny
      dn = δ[i,j]
      if dn > 0
        px[s:s+dn-1,1] .+= config.grid.x[i,j]
        px[s:s+dn-1,2] .+= config.grid.y[i,j]
        s += dn
      end
    end
  end
  remove_particles!(species, config.grid.Δh, (i,j) -> begin
    if δ[i,j] < 0
      δ[i,j] += 1
      return true
    end
    return false
  end)
  species.np += n
end
