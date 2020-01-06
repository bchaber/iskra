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
Δt = 50ns       # time step
Lx = nx*Δh      # domain length in x direction
Ly = ny*Δh      # domain length in y direction
############################################
xs, ys = 0m:Δh:Lx, 0m:Δh:Ly
sx, sv = [0 Lx; 0 Ly], [0 0; 0 0]
O  = create_fluid_species( "O", 0.0, 0qe, 8mp, nx+1, ny+1)
iO = create_fluid_species("O+", 0.0,+1qe, 8mp, nx+1, ny+1)
e  = create_fluid_species("e-", 1.0,-1qe, 1me, nx+1, ny+1)
using Chemistry
import RegularGrid, FiniteDifferenceMethod, ParticleInCell
config.grid    = RegularGrid.create_uniform_grid(xs, ys)
config.cells   = RegularGrid.create_staggered_grid(config.grid)
config.solver  = FiniteDifferenceMethod.create_poisson_solver(config.grid)
config.pusher  = ParticleInCell.create_boris_pusher()
config.species = [O, iO, e]

#σ = CrossSection(0:0.3:1.5, [0, 0.1e-4, 0.4e-4, 0.5e-4, 0.7e-4, 0.9e-4])
σ(E) = .1
chemistry = chemical(@reactions begin
    σ, e + O --> 2e + iO
end)
config.interactions = [chemistry]

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
  Diagnostics.open_container("03-field")
  Diagnostics.open_container("03-particle")
end

function ParticleInCell.after_loop(it)
  Diagnostics.save_diagnostic("E",   "03-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("ϕ",   "03-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("nO",  "03-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("nO+", "03-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("ne-", "03-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("ΔnO", "03-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("ΔnO+","03-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("Δne-","03-field",   it, Δt*it-Δt)
end

function ParticleInCell.exit_loop()
  Diagnostics.close_container("03-field")
  Diagnostics.close_container("03-particle")
end
ParticleInCell.init(ParticleInCell.DensitySource(1e6δ, config.grid), O, Δt)
ParticleInCell.init(ParticleInCell.DensitySource(1e4δ, config.grid), e, Δt)
@time ParticleInCell.solve(config, Δt, ts, ε0)