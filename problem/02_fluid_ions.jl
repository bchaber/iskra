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
Δt = 5000ns     # time step
Lx = nx*Δh      # domain length in x direction
Ly = ny*Δh      # domain length in y direction
############################################
xs, ys = 0m:Δh:Lx, 0m:Δh:Ly
sx, sv = [0 Lx; 0 Ly], [0 0; 0 0]
iO  = create_fluid_species("O+", 1.0, +1qe, 8mp, nx+1, ny+1)
import RegularGrid, FiniteDifferenceMethod, ParticleInCell
config.grid    = RegularGrid.create_uniform_grid(xs, ys)
config.cells   = RegularGrid.create_staggered_grid(config.grid)
config.solver  = FiniteDifferenceMethod.create_poisson_solver(config.grid, ε0)
config.species = [iO]
############################################
nx, ny = size(config.grid)
mx, my = size(config.cells)
xx, yy = config.grid.coords
δ = @. exp(-(xx-0.5Lx)^2/0.03Lx -(yy-0.5Ly)^2/0.03Ly)
εr  = ones(mx, my, 1)
bcs = zeros(Int8, nx, ny, 1)
bcs[ nx,  1, 1] = 1
bcs[ nx, ny, 1] = 2
set_permittivity(εr)
create_electrode(bcs .== 1, config; σ=30ε0)
create_electrode(bcs .== 2, config; fixed=true)
############################################
import ParticleInCell
using Diagnostics
using XDMF

function ParticleInCell.after_loop(i, t, dt)
  cd("/tmp")
  new_iteration("02_fluid_ions", i, t, dt) do it
    save_record(it, "nO+")
    save_record(it, "dO+")
    save_record(it, "vO+")
    save_record(it, "phi")
    save_record(it, "E")
  end
end

function ParticleInCell.exit_loop()
  println("Exporting to XDMF...")
  cd("/tmp/02_fluid_ions")
  xdmf(1:ts) do it
  	save_fields(it,  "xdmf/fields.xdmf")
  end
end

ParticleInCell.init(ParticleInCell.DensitySource(1e5δ, config.grid), iO, Δt)
@time ParticleInCell.solve(config, Δt, ts)