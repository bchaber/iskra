using RegularGrid

mutable struct DensitySource{D}
  δ :: Array{Float64, D}
  grid :: CartesianGrid{D}
end

mutable struct MaxwellianSource{D, V}
   rate :: Float64
   x :: AbstractArray{Float64, 2}
   v :: AbstractArray{Float64, 2}
end

function sample!(config :: MaxwellianSource{D, V}, species :: KineticSpecies{D, V}, Δt) where {D, V}
  np = species.np
  cx, cv = config.x, config.v
  px, pv = @views species.x[1+np:end,:], species.v[1+np:end,:]

  n = minimum([size(px, 1), floor(Integer, config.rate*Δt)])
  Dx = reshape(cx[:, 1], 1, :)
  dx = reshape(cx[:, 2], 1, :)
  Dv = reshape(cv[:, 1], 1, :)
  dv = reshape(cv[:, 2], 1, :)
  px[1:n, :] =  rand(n,D) .* repeat(dx, n) .+ repeat(Dx, n)
  pv[1:n, :] = randn(n,V) .* repeat(dv, n) .+ repeat(Dv, n)
  species.np += n
end

function sample!(config :: DensitySource{D}, species :: FluidSpecies{D}, Δt) where D
  species.n .+= config.δ
end