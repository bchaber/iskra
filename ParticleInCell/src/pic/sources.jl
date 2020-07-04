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