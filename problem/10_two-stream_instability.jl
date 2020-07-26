import Random
Random.seed!(0)
############################################
include("configuration.jl")
config = Config()
############################################
include("units_and_constants.jl")
############################################
nx = 20         # number of nodes in x direction
ny = 2          # number of nodes in y direction
ts = 200        # number of time steps
Δh = 5cm        # cell size
Δt = 20ns       # time step
Lx = nx*Δh      # domain length in x direction
Ly = ny*Δh      # domain length in y direction
############################################
xs, ys = 0m:Δh:Lx, 0:Δh:Ly
sx, sv = [0 Lx; 0 Ly], [0 -0.05Δh/Δt; 0 0]
sv = [0 -0.05Δh/Δt]
e = create_kinetic_species("e-", 20_000,-1qe, 1me, 1);
γ = create_gamma_ionization_source(1/Δt, sx, sv);

using RegularGrid, FiniteDifferenceMethod, ParticleInCell
config.grid    = create_uniform_grid(xs)
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
using Diagnostics
using XDMF

function ParticleInCell.enter_loop(config)
end

function ParticleInCell.after_loop(i, t, dt)
  cd("/tmp")
  new_iteration("10_two-stream_instability", i, t, dt) do it
    save_records(it, "e-/")
    save_record(it, "rho")
    save_record(it, "phi")
  end
end

function ParticleInCell.exit_loop()
  println("Exporting to XDMF...")
  cd("/tmp/10_two-stream_instability")
  xdmf(1:ts) do it
  	save_species(it, "xdmf/electrons.xdmf", "e-")
  	save_fields(it,  "xdmf/fields.xdmf")
  end
end

init(γ, e, Δt)
@time solve(config, Δt, ts)