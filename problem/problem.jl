import Random
Random.seed!(0)
############################################
include("configuration.jl")
config = Config()
############################################
include("units_and_constants.jl")
############################################
ϕ0 =-1V         # reference potential
ϕp =+1V         # wall potential
Ti = 0.1eV      # ion temperature in eV
v_drift = 5mps # ion injection velocity
v_therm = 0#sqrt(2qe*Ti/100u)  # thermal velocity with Ti in eV
nx = 20         # number of nodes in x direction
ny = 20         # number of nodes in y direction
ts = 200        # number of time steps
Δh =  5cm       # cell size
Δt = 100ns      # time step
nn = nx*ny      # total number of nodes
Lx = nx*Δh      # domain length in x direction
Ly = ny*Δh      # domain length in y direction
############################################
xs, ys = 0m:Δh:Lx, 0m:Δh:Ly
sx, sv = [0 Lx; 0 Ly], [v_therm v_drift; v_therm 0]
e  = create_kinetic_species("e-", 20_000,-1qe, 1me, 10_0000)
iO = create_kinetic_species("O+", 20_000,+1qe,  1u,  1_0000)
#chem = create_chemical_reactions("O")
γe = create_gamma_ionization_source( e, .5/Δt, sx, sv)

import RegularGrid, FiniteDifferenceMethod, ParticleInCell
config.grid    = RegularGrid.create_uniform_grid(xs, ys)
config.cells   = RegularGrid.create_staggered_grid(config.grid)
config.solver  = FiniteDifferenceMethod.create_poisson_solver(config.grid)
config.pusher  = ParticleInCell.create_boris_pusher()
config.species = [iO, e]
config.sources = [γe]

sO = ParticleInCell.create_maxwellian_source(iO, 10/Δt, sx, sv)
############################################
nx, ny = size(config.grid)
mx, my = size(config.cells)
εr  = ones(mx, my, 1)
bcs = zeros(Int8, nx, ny)
inbox = (0.4m .<= config.cells.x .<= 0.6m) .&
        (-.5m .<= config.cells.y .<= 1.5m)
inbox = reshape(inbox, mx, my, 1)
bcs[1, 1:ny] .= 1
bcs[nx,1:ny] .= 2
εr[inbox] .= 1.
config.cells["εr"] = εr
set_permittivity(εr)
add_electrode(bcs .== 1, ϕ0)
add_electrode(bcs .== 2, ϕp)
############################################
import ParticleInCell
import Diagnostics

function ParticleInCell.enter_loop()
  Diagnostics.open_container("problem-field")
  Diagnostics.open_container("problem-particle")
  Diagnostics.open_container("problem-cells")
end

function ParticleInCell.after_loop(it)
  Diagnostics.save_diagnostic("εr",  "problem-cells",   it)
  Diagnostics.save_diagnostic("ρ",   "problem-field",   it)
  Diagnostics.save_diagnostic("E",   "problem-field",   it)
  Diagnostics.save_diagnostic("ϕ",   "problem-field",   it)
  Diagnostics.save_diagnostic("nO+", "problem-field",   it)
  Diagnostics.save_diagnostic("ne-", "problem-field",   it)
  Diagnostics.save_diagnostic("pvO+","problem-particle",it)
  Diagnostics.save_diagnostic("pve-","problem-particle",it)
end

function ParticleInCell.exit_loop()
  Diagnostics.close_container("problem-field")
  Diagnostics.close_container("problem-particle")
  Diagnostics.close_container("problem-cells")
end

ParticleInCell.init(sO, Δt)
@time ParticleInCell.solve(config, Δt, ts, ε0)
