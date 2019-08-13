module Source
export create_particles!

include("pic/macroparticle.jl")

mutable struct MaxwellianSourceConfiguration
   species
   n :: Integer
   x :: AbstractArray{Float64,2}
   v :: AbstractArray{Float64,2}
end

function create_maxwellian_source(sps, n, x, v)
  MaxwellianSourceConfiguration(sps, n, x, v)
end

function create_particles!(config :: MaxwellianSourceConfiguration, it)
  part = config.species
  cx, cv = config.x, config.v
  px, pv = @views part.x[1+part.np:end,:], part.v[1+part.np:end,:]

  if it > 1
    return 0
  end

  n = minimum([config.n,size(px,1)])
  # insert particles randomly distributed in y and in the first cell
  px[1:n,1]=randn(n,1)*cx[1,2]/10 .+ cx[1,2]/2 # x position
  px[1:n,2]=randn(n,1)*cx[2,2]/10 .+ cx[2,2]/2 # y position
  # sample Maxwellian in x and y, add drift velocity in x
  pv[1:n,1]=1.0*(-1.5.+rand(n,1).+rand(n,1).+rand(n,1))*cv[1,1] .+ cv[1,2]
  pv[1:n,2]=0.5*(-1.5.+rand(n,1).+rand(n,1).+rand(n,1))*cv[2,1] .+ cv[2,2]
  part.np += n

  return n
end
end
