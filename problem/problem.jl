import Random
Random.seed!(0)
############################################
include("configuration.jl")
config = Config()
############################################
include("units_and_constants.jl")
############################################
n0 = 1e12       # density in #/m^3
ϕ0 =-1V         # reference potential
ϕp =+1V         # wall potential
Te = 1.0eV      # electron temperature in eV
Ti = 0.1eV      # ion temperature in eV
v_drift = 0#-1mps # ion injection velocity
v_th = 0#sqrt(2qe*Ti/100u)  # thermal velocity with Ti in eV
λD = sqrt(ɛ0*Te/(n0*qe)) # Debye length
# set simulation domain
nx = 50         # number of nodes in x direction
ny = 50         # number of nodes in y direction
ts = 75         # number of time steps
Δh = 1cm        # cell size
Δt = 30μs       # time step
#Δt = 100ns     # time step
nn = nx*ny      # total number of nodes
Lx = nx*Δh      # domain length in x direction
Ly = ny*Δh      # domain length in y direction
############################################
xs, ys = 0m:Δh:Lx, 0m:Δh:Ly
sx, sv = [0 Lx; 0 Ly], [v_th v_drift; v_th 0]
e  = create_species("e-", 20_000,-1qe, 1me, 100)
iO = create_species("O+", 20_000,+1qe,100u, 100)

import RegularGrid, FiniteDifferenceMethod
import ParticleInCell: Source, Pusher
config.grid    = RegularGrid.create_uniform_grid(xs, ys)
config.solver  = FiniteDifferenceMethod.create_poisson_solver(config.grid)
config.source  = Source.create_maxwellian_source(iO, sx, sv)
config.pusher  = Pusher.create_boris_pusher()
config.species = [iO]
############################################
nx, ny = size(config.grid)
bcs = zeros(Int8, nx, ny)
inbox = (0.02m .<= config.grid.x .<= 0.05m) .&
        (0.02m .<= config.grid.y .<= 0.04m)
bcs[inbox].= 5;
bcs[ 1, 1:5] .= 1
bcs[nx, 1:5] .= 2

add_electrode(bcs .== 1,-1V)
add_electrode(bcs .== 2,+1V)
############################################
import ParticleInCell
import Diagnostics

function ParticleInCell.enter_loop()
  Diagnostics.open_container("problem-field")
  Diagnostics.open_container("problem-particle")
end

function ParticleInCell.after_loop(it)
  Diagnostics.save_diagnostic("ρ","problem-field",it)
  Diagnostics.save_diagnostic("ϕ","problem-field",it)
  Diagnostics.save_diagnostic("E","problem-field",it)
  Diagnostics.save_diagnostic("pvO+","problem-particle",it)
  Diagnostics.save_diagnostic("pEO+","problem-particle",it)
end

function ParticleInCell.exit_loop()
  Diagnostics.close_container("problem-field")
  Diagnostics.close_container("problem-particle")
end

@time ParticleInCell.solve(config, Δt, ts)
