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
using Diagnostics
using XDMF

function ParticleInCell.after_loop(i, t, dt)
  new_iteration("07_boundaries", i, t, dt) do it
    save_records(it, "e-/")
    save_record(it, "ne-")
    save_record(it, "phi")
    save_record(it, "E")
  end
end

function ParticleInCell.exit_loop()
  println("Exporting to XDMF...")
  xdmf("07_boundaries", 1:ts, "electrons.xdmf") do it, f, t
  	save_species(["e-"], it, f, t)
  end
  xdmf("07_boundaries", 1:ts, "fields.xdmf") do it, f, t
  	save_fields(it, f, t)
  end
end

ParticleInCell.init(γ, e, Δt)
@time ParticleInCell.solve(config, Δt, ts)