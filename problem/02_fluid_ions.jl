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
config.solver  = FiniteDifferenceMethod.create_poisson_solver(config.grid)
config.species = [iO]
############################################
nx, ny = size(config.grid)
mx, my = size(config.cells)
xx, yy = config.grid.x, config.grid.y
δ = @. exp(-(xx-0.5Lx)^2/0.03Lx -(yy-0.5Ly)^2/0.03Ly)
εr  = ones(mx, my, 1)
bcs = zeros(Int8, nx, ny, 1)
bcs[ nx,  1, 1] = 1
bcs[ nx, ny, 1] = 2
set_permittivity(εr)
create_electrode(bcs .== 1, config.solver, config.grid; fixed=true, V=+1V)
create_electrode(bcs .== 2, config.solver, config.grid; fixed=true, V=-1V)
############################################
import ParticleInCell
import Diagnostics

function ParticleInCell.enter_loop()
  Diagnostics.open_container("02-field")
end

function ParticleInCell.after_loop(it)
  Diagnostics.save_diagnostic("E",   "02-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("ϕ",   "02-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("nO+", "02-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("vO+", "02-field",   it, Δt*it-Δt)
end

function ParticleInCell.exit_loop()
  Diagnostics.close_container("02-field")
end

ParticleInCell.init(ParticleInCell.DensitySource(1e5δ, config.grid), iO, Δt)
@time ParticleInCell.solve(config, Δt, ts, ε0)