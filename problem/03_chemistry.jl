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
iO = create_kinetic_species("O+", 300_000,+1qe, 8mp, 100)
e  = create_kinetic_species("e-", 300_000,-1qe, 1me, 100)
using Chemistry
import RegularGrid, FiniteDifferenceMethod, ParticleInCell
config.grid    = RegularGrid.create_uniform_grid(xs, ys)
config.cells   = RegularGrid.create_staggered_grid(config.grid)
config.solver  = FiniteDifferenceMethod.create_poisson_solver(config.grid)
config.pusher  = ParticleInCell.create_boris_pusher()
config.species = [O, iO, e]

#σ = CrossSection(0:0.3:1.5, [0, 0.1e-4, 0.4e-4, 0.5e-4, 0.7e-4, 0.9e-4])
σ(E) = 100
config.chemistry = @reactions begin
    σ, e + O --> 2e + iO
end
println(config.chemistry)

############################################
nx, ny = size(config.grid)
mx, my = size(config.cells)
xx, yy = config.grid.x, config.grid.y
δ = @. exp(-(xx-0.5Lx)^2/0.03Lx -(yy-0.5Ly)^2/0.03Ly)
εr  = ones(mx, my, 1)
bcs = zeros(Int8, nx, ny)
bcs[ nx,  1] = 1
bcs[ nx, ny] = 2
set_permittivity(εr)
add_electrode(bcs .== 1, +1V)
add_electrode(bcs .== 2, -1V)
############################################
import ParticleInCell
import Diagnostics

function ParticleInCell.enter_loop()
  Diagnostics.open_container("problem-field")
  Diagnostics.open_container("problem-particle")
end

function ParticleInCell.after_loop(it)
  Diagnostics.save_diagnostic("E",   "problem-field",   it)
  Diagnostics.save_diagnostic("ϕ",   "problem-field",   it)
  Diagnostics.save_diagnostic("nO",  "problem-field",   it)
  Diagnostics.save_diagnostic("nO+", "problem-field",   it)
  Diagnostics.save_diagnostic("ne-", "problem-field",   it)
  Diagnostics.save_diagnostic("ΔnO", "problem-field",   it)
  Diagnostics.save_diagnostic("ΔnO+","problem-field",   it)
  Diagnostics.save_diagnostic("Δne-","problem-field",   it)
  Diagnostics.save_diagnostic("pvO+","problem-particle",it)
  Diagnostics.save_diagnostic("pve-","problem-particle",it)
end

function ParticleInCell.exit_loop()
  Diagnostics.close_container("problem-field")
  Diagnostics.close_container("problem-particle")
end
ParticleInCell.init(ParticleInCell.DensitySource(1e5δ, config.grid), O, Δt)
ParticleInCell.init(ParticleInCell.DensitySource(1e5δ, config.grid), e, Δt)
@time ParticleInCell.solve(config, Δt, ts, ε0)