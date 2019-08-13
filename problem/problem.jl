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
v_drift =  5mps # ion injection velocity
v_th = 0#sqrt(2qe*Ti/100u)  # thermal velocity with Ti in eV
λD = sqrt(ɛ0*Te/(n0*qe)) # Debye length
# set simulation domain
nx = 20         # number of nodes in x direction
ny = 20         # number of nodes in y direction
ts = 75         # number of time steps
Δh =  5cm       # cell size
#Δt = 30μs       # time step
Δt = 100ns     # time step
nn = nx*ny      # total number of nodes
Lx = nx*Δh      # domain length in x direction
Ly = ny*Δh      # domain length in y direction
############################################
xs, ys = 0m:Δh:Lx, 0m:Δh:Ly
sx, sv = [0 Lx; 0 Ly], [v_th v_drift; v_th 0]
e  = create_species("e-", 20_000,-1qe, 1me, 100)
iO = create_species("O+", 20_000,+1qe,  1u, 100)
#chem = create_chemical_reactions("O")
γe = create_gamma_ionization_source( e, 500, sx, sv)
γO = create_gamma_ionization_source(iO,  10, sx, sv)

import RegularGrid, FiniteDifferenceMethod
import ParticleInCell: Source, Pusher
config.grid    = RegularGrid.create_uniform_grid(xs, ys)
config.solver  = FiniteDifferenceMethod.create_poisson_solver(config.grid)
config.pusher  = Pusher.create_boris_pusher()
config.species = [iO, e]
config.sources = [γe, γO]
############################################
nx, ny = size(config.grid)
bcs = zeros(Int8, nx, ny)
inbox = (0.02m .<= config.grid.x .<= 0.05m) .&
        (0.02m .<= config.grid.y .<= 0.04m)
bcs[inbox].= 5;
bcs[1,  1] = 1
bcs[1, ny] = 2

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
  Diagnostics.save_diagnostic("ρ",   "problem-field",   it)
  Diagnostics.save_diagnostic("nO+", "problem-field",   it)
  Diagnostics.save_diagnostic("ne-", "problem-field",   it)
  Diagnostics.save_diagnostic("pvO+","problem-particle",it)
  Diagnostics.save_diagnostic("pve-","problem-particle",it)
end

function ParticleInCell.exit_loop()
  Diagnostics.close_container("problem-field")
  Diagnostics.close_container("problem-particle")
end

@time ParticleInCell.solve(config, Δt, ts)
