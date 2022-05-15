import Random

mutable struct MaxwellianSource{D, V}
   rate :: Float64
   wx :: SVector{D, Float64}
   dx :: SVector{D, Float64}
   wv :: SVector{V, Float64}
   dv :: SVector{V, Float64}
   MaxwellianSource{D,V}(rate::Float64, wx::Vector{Float64}, wv::Vector{Float64};
    dx=nothing, dv=nothing) where {D,V} =
     new{D,V}(rate,
              wx, isnothing(dx) ? zeros(D) : dx,
              wv, isnothing(dv) ? zeros(V) : dv)
end

function sample!(config :: MaxwellianSource{D, V}, species :: KineticSpecies{D, V}, Δt) where {D, V}
  #println("Sampling with the following seed: ", Random.GLOBAL_SEED)
  wx, wv = config.wx, config.wv
  dx, dv = config.dx, config.dv
  px, pv = species.x, species.v
  n = min(length(species.id), floor(Integer, config.rate*Δt))
  for p in (species.np + 1) : (species.np + n)
    px[p] =  rand(D) .* wx .+ dx
    pv[p] = randn(V) .* wv .+ dv
  end
  species.np += n
end