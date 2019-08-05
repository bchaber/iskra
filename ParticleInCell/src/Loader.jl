module Loader

include("pic/macroparticle.jl")

function setup(mnp, xe, ye, vth, vdr)
  global xmax
  global ymax
  global vthermal
  global vdrift
  global npmax
  xmax, ymax = xe, ye
  vthermal = vth
  vdrift = vdr
  npmax = mnp
end

function load_particles(species, iteration, np_limit)
  q, m, np2c = species.q, species.m, species.np2c
  name = species.name
  if iteration > 1
    return Species(name, 0)
  end

  n = minimum([np_limit npmax])
  x = zeros(n, 2)
  v = zeros(n, 2)
  # insert particles randomly distributed in y and in the first cell
  x[:,1]=randn(n,1)*xmax/10 .+ xmax/2 # x position
  x[:,2]=randn(n,1)*ymax/10 .+ ymax/2 # y position
  # sample Maxwellian in x and y, add drift velocity in x
  v[:,1]=1.0*(-1.5.+rand(n,1).+rand(n,1).+rand(n,1))*vthermal .- vdrift
  v[:,2]=0.5*(-1.5.+rand(n,1).+rand(n,1).+rand(n,1))*vthermal

  return Species(x, v, q, m, name, np2c, n)
end
end
