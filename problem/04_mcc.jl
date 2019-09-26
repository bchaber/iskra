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
Δt = .5ns        # time step
Lx = nx*Δh      # domain length in x direction
Ly = ny*Δh      # domain length in y direction
############################################
xs, ys = 0m:Δh:Lx, 0m:Δh:Ly
sx, sv = [0 Lx; 0 Ly], [0 0; 0 0]
O  = create_fluid_species("O", 1.0, 0qe, 8mp, nx+1, ny+1)
e  = create_kinetic_species("e-", 20_000,-1qe, 1me, 1_000)
iO = create_kinetic_species("O+", 20_000,+1qe, 8mp, 1_000)
using Chemistry
import RegularGrid, FiniteDifferenceMethod, ParticleInCell
config.grid    = RegularGrid.create_uniform_grid(xs, ys)
config.cells   = RegularGrid.create_staggered_grid(config.grid)
config.solver  = FiniteDifferenceMethod.create_poisson_solver(config.grid)
config.pusher  = ParticleInCell.create_boris_pusher()
config.species = [e, O, iO]

#σ = CrossSection(0:0.3:1.5, [0, 0.1e-7, 0.4e-7, 0.5e-7, 0.7e-7, 0.9e-7])
σ(g) = 5e-16g
collisions = ParticleInCell.mcc(@reactions begin
    #σ, e + O --> O + e
    σ, e + O --> iO + 2e
end)
config.interactions = [collisions]
############################################
nx, ny = size(config.grid)
mx, my = size(config.cells)
xx, yy = config.grid.x, config.grid.y
δ = @. exp(-(xx-0.5Lx)^2/0.02Lx -(yy-0.5Ly)^2/0.02Ly)
εr  = ones(mx, my, 1)
bcs = zeros(Int8, nx, ny)
bcs[ nx,  1] = 1
bcs[ nx, ny] = 2
set_permittivity(εr)
add_electrode(bcs .== 1, +1e3V)
add_electrode(bcs .== 2, -1e3V)
############################################
import ParticleInCell
import Diagnostics

function ParticleInCell.enter_loop()
  Diagnostics.open_container("04-field")
  Diagnostics.open_container("04-particle")
end

function ParticleInCell.after_loop(it)
  Diagnostics.save_diagnostic("E",   "04-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("ϕ",   "04-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("nO",  "04-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("ne-", "04-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("nO+", "04-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("pve-","04-particle",it, Δt*it-Δt)
  Diagnostics.save_diagnostic("pvO+","04-particle",it, Δt*it-Δt)
end

function ParticleInCell.exit_loop()
  Diagnostics.close_container("04-field")
  Diagnostics.close_container("04-particle")
end
ParticleInCell.init(ParticleInCell.DensitySource(5e7δ, config.grid), e, Δt)
ParticleInCell.init(ParticleInCell.DensitySource(5e7δ, config.grid), O, Δt)
@time ParticleInCell.solve(config, Δt, 66, ε0)