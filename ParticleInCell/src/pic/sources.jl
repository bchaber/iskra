using RegularGrids

mutable struct DensitySource{D}
  δ :: Array{Float64, D}
  grid :: CartesianGrid{D}
end

mutable struct MaxwellianSource{D, V}
   rate :: Float64
   wx :: Matrix{Float64}
   dx :: Matrix{Float64}
   wv :: Matrix{Float64}
   dv :: Matrix{Float64}
   MaxwellianSource{1,V}(rate::Float64, wx::Float64, wv::Matrix{Float64}; dx=zero(Float64), dv=nothing) where {V} =
    MaxwellianSource{1,V}(rate, [wx], wv; dx=[dx], dv=dv)
   MaxwellianSource{1,V}(rate::Float64, wx::Vector{Float64}, wv::Matrix{Float64}; dx=zeros(1),dv=nothing) where {V} =
    MaxwellianSource{1,V}(rate, wx[:,:], wv; dx=dx[:,:], dv=dv)
   MaxwellianSource{D,V}(rate::Float64, wx::Matrix{Float64}, wv::Matrix{Float64}; dx=nothing, dv=nothing) where {D,V} =
     new{D,V}(rate,
              wx, isnothing(dx) ? zeros(1, D) : dx,
              wv, isnothing(dv) ? zeros(1, V) : dv)
end

function sample!(config :: MaxwellianSource{D, V}, species :: KineticSpecies{D, V}, Δt) where {D, V}
  np = species.np
  wx, wv = config.wx, config.wv
  dx, dv = config.dx, config.dv
  px, pv = @views species.x[1+np:end,:], species.v[1+np:end,:]

  n = minimum([size(px, 1), floor(Integer, config.rate*Δt)])
  px[1:n, :] =  rand(n,D) .* repeat(wx, n) .+ repeat(dx, n)
  pv[1:n, :] = randn(n,V) .* repeat(wv, n) .+ repeat(dv, n)
  species.np += n
end

function sample!(config :: DensitySource{D}, species :: FluidSpecies{D}, Δt) where D
  species.n .+= config.δ
end