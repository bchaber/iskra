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

using RegularGrid, FiniteDifferenceMethod, ParticleInCell
config.grid    = create_uniform_grid(xs, ys)
config.cells   = create_staggered_grid(config.grid)
config.solver  = create_poisson_solver(config.grid, ε0)
config.pusher  = create_boris_pusher()
config.species = [e]
############################################
nx, ny = size(config.grid)
mx, my = size(config.cells)
εr  = ones(mx, my, 1)
bcs = zeros(Int8, nx, ny, 1)
bcs[ nx, ny, 1] = 1
bcs[ nx,  1, 1] = 2
set_permittivity(εr)
create_electrode(bcs .== 1, config; σ=-1ε0)
create_electrode(bcs .== 2, config; fixed=true)
############################################
import ParticleInCell
import Diagnostics

function ParticleInCell.after_loop(i, t, dt)
  Diagnostics.new_iteration("01_single_particle", i, t, dt) do it
    Diagnostics.save_diagnostic(it, "e-/positionOffset/x")
    Diagnostics.save_diagnostic(it, "e-/positionOffset/y")
    Diagnostics.save_diagnostic(it, "e-/positionOffset/z")
    Diagnostics.save_diagnostic(it, "e-/position")
    Diagnostics.save_diagnostic(it, "e-/momentum")
    Diagnostics.save_diagnostic(it, "e-/weighting")
    Diagnostics.save_diagnostic(it, "e-/charge")
    Diagnostics.save_diagnostic(it, "e-/mass")
    Diagnostics.save_diagnostic(it, "e-/id")
    Diagnostics.save_diagnostic(it, "rho")
    Diagnostics.save_diagnostic(it, "phi")
  end
end

init(γ, e, Δt)
@time solve(config, Δt, ts)