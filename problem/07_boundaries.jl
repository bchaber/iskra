import Random
Random.seed!(0)
############################################
include("configuration.jl")
config = Config()
############################################
include("units_and_constants.jl")
############################################
nx = 10         # number of nodes in x direction
ny = 10         # number of nodes in y direction
ts = 200        # number of time steps
Δh = 10cm       # cell size
Δt = 10ns       # time step
Lx = nx*Δh      # domain length in x direction
Ly = ny*Δh      # domain length in y direction
############################################
xs, ys = 0m:Δh:Lx, 0m:Δh:Ly
sx, sv = [Lx/4 Lx/4; Ly/4 Ly/2], [0Δh/Δt 0Δh/Δt; 0Δh/Δt 0Δh/Δt]
e = create_kinetic_species("e-", 50_000,-1qe, 1me, 50e3)
γ = create_gamma_ionization_source(20_000/Δt, sx, sv)

import RegularGrid, FiniteDifferenceMethod, ParticleInCell
config.grid    = RegularGrid.create_uniform_grid(xs, ys)
config.cells   = RegularGrid.create_staggered_grid(config.grid)
config.solver  = FiniteDifferenceMethod.create_poisson_solver(config.grid, ε0)
config.pusher  = ParticleInCell.create_boris_pusher()
config.species = [e]
############################################
nx, ny = size(config.grid)
bcs = zeros(Int8, nx, ny, 1)
bcs[ 1, 2:ny-1, 1] .= 1
bcs[nx, 5:7,    1] .= 2
bcs[nx-1,  1,   1]  = 3
bcs[nx-1, ny,   1]  = 3
bcs[6:8,  5:7,  1] .= 4
driven   = create_electrode(bcs .== 1, config.solver, config.grid; σ=1ε0)
floating = create_electrode(bcs .== 2, config.solver, config.grid)
grounded = create_electrode(bcs .== 3, config.solver, config.grid; fixed=true)
reflecting = ParticleInCell.create_reflective_surface()
config.tracker = ParticleInCell.create_surface_tracker(bcs,
	[driven, grounded, floating, reflecting], Δh, Δt)
############################################
import ParticleInCell
import Diagnostics

function ParticleInCell.enter_loop()
  Diagnostics.open_container("07-field")
  Diagnostics.open_container("07-particle")
end

function ParticleInCell.after_loop(it)
  Diagnostics.save_diagnostic("E",   "07-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("ϕ",   "07-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("ne-", "07-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("pve-","07-particle",it, Δt*it-Δt)
end

function ParticleInCell.exit_loop()
  Diagnostics.close_container("07-field")
  Diagnostics.close_container("07-particle")
end

ParticleInCell.init(γ, e, Δt)
@time ParticleInCell.solve(config, Δt, ts)