module Source
export create_particles

include("pic/macroparticle.jl")

mutable struct MaxwellianSourceConfiguration
   sp
   x :: AbstractArray{Float64,2}
   v :: AbstractArray{Float64,2}
end

function create_maxwellian_source(sp, x, v)
  MaxwellianSourceConfiguration(sp, x, v)
end

function create_particles(config :: MaxwellianSourceConfiguration, it, n)
  q, m, np2c = config.sp.q, config.sp.m, config.sp.np2c
  cx, cv = config.x, config.v
  name = config.sp.name

  if it > 1
    return Species(name, 0)
  end

  x = zeros(n, 2)
  v = zeros(n, 2)
  # insert particles randomly distributed in y and in the first cell
  x[:,1]=randn(n,1)*cx[1,2]/10 .+ cx[1,2]/2 # x position
  x[:,2]=randn(n,1)*cx[2,2]/10 .+ cx[2,2]/2 # y position
  # sample Maxwellian in x and y, add drift velocity in x
  v[:,1]=1.0*(-1.5.+rand(n,1).+rand(n,1).+rand(n,1))*cv[1,1] .+ cv[1,2]
  v[:,2]=0.5*(-1.5.+rand(n,1).+rand(n,1).+rand(n,1))*cv[2,1] .+ cv[2,2]

  return Species(x, v, q, m, name, np2c, n)
end
end
