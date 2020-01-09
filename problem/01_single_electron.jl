import Random
Random.seed!(0)
############################################
include("configuration.jl")
config = Config()
############################################
include("units_and_constants.jl")
############################################
nx = 20         # number of nodes in x direction
ny = 20         # number of nodes in y direction
ts = 200        # number of time steps
Δh = 5cm        # cell size
Δt = 20ns       # time step
Lx = nx*Δh      # domain length in x direction
Ly = ny*Δh      # domain length in y direction
############################################
xs, ys = 0m:Δh:Lx, 0m:Δh:Ly
sx, sv = [0 Lx; 0 Ly], [0 -.1Δh/Δt; 0 0]
e = create_kinetic_species("e-", 20_000,-1qe, 1me, 1)
γ = create_gamma_ionization_source(1/Δt, sx, sv)

import RegularGrid, FiniteDifferenceMethod, ParticleInCell
config.grid    = RegularGrid.create_uniform_grid(xs, ys)
config.cells   = RegularGrid.create_staggered_grid(config.grid)
config.solver  = FiniteDifferenceMethod.create_poisson_solver(config.grid)
config.pusher  = ParticleInCell.create_boris_pusher()
config.species = [e]
############################################
nx, ny = size(config.grid)
mx, my = size(config.cells)
εr  = ones(mx, my, 1)
bcs = zeros(Int8, nx, ny, 1)
bcs[ nx,  1, 1] = 1
bcs[ nx, ny, 1] = 2
set_permittivity(εr)
create_electrode(bcs .== 1, config.solver, config.grid; fixed=true, ϕ=+1V)
create_electrode(bcs .== 2, config.solver, config.grid; fixed=true, ϕ=-1V)
############################################
import ParticleInCell
import Diagnostics

function ParticleInCell.enter_loop()
  Diagnostics.open_container("01-field")
  Diagnostics.open_container("01-particle")
end

function ParticleInCell.after_loop(it)
  Diagnostics.save_diagnostic("E",   "01-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("ϕ",   "01-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("ne-", "01-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("pve-","01-particle",it, Δt*it-Δt)
end

function ParticleInCell.exit_loop()
  Diagnostics.close_container("01-field")
  Diagnostics.close_container("01-particle")
end

ParticleInCell.init(γ, e, Δt)
@time ParticleInCell.solve(config, Δt, ts, ε0)