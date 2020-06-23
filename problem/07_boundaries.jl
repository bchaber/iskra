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
reflecting = ParticleInCell.create_reflective_surface()
driven   = create_electrode(bcs .== 1, config; σ=1ε0)
floating = create_electrode(bcs .== 2, config)
grounded = create_electrode(bcs .== 3, config; fixed=true)
ParticleInCell.track_surface!(config.tracker, bcs .== 4, reflecting)
############################################
import ParticleInCell
import Diagnostics

function ParticleInCell.after_loop(i, t, dt)
  Diagnostics.new_iteration("07_boundaries", i, t, dt) do it
    Diagnostics.save_diagnostic(it, "ne-")
    Diagnostics.save_diagnostic(it, "e-/positionOffset/x")
    Diagnostics.save_diagnostic(it, "e-/positionOffset/y")
    Diagnostics.save_diagnostic(it, "e-/positionOffset/z")
    Diagnostics.save_diagnostic(it, "e-/position")
    Diagnostics.save_diagnostic(it, "e-/momentum")
    Diagnostics.save_diagnostic(it, "e-/weighting")
    Diagnostics.save_diagnostic(it, "e-/charge")
    Diagnostics.save_diagnostic(it, "e-/mass")
    Diagnostics.save_diagnostic(it, "e-/id")
    Diagnostics.save_diagnostic(it, "phi")
    Diagnostics.save_diagnostic(it, "E")
  end
end

ParticleInCell.init(γ, e, Δt)
@time ParticleInCell.solve(config, Δt, ts)