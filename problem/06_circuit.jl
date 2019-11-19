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
e  = create_kinetic_species("e-", 20_000,-1qe, 1me, 1)
iO = create_kinetic_species("O+", 20_000,+1qe, 8mp, 1)
using Chemistry, Circuit
import RegularGrid, FiniteDifferenceMethod, ParticleInCell
config.grid    = RegularGrid.create_uniform_grid(xs, ys)
config.cells   = RegularGrid.create_staggered_grid(config.grid)
config.solver  = FiniteDifferenceMethod.create_poisson_solver(config.grid)
config.pusher  = ParticleInCell.create_boris_pusher()
config.species = [e, O, iO]
cir = rlc(@netlist begin
    R1, R1₊, R1₋, 50
    C,  R1₋, GND, 200,
    V,  R1₊, GND, sin(2π*f*t)
end)
σ = CrossSection(1.5e7:7.5e7:9e7, [0.1, 1.0])
collisions = mcc(@reactions begin
    σ, e + O --> O + e
end)
config.interactions = [collisions]
############################################
nx, ny = size(config.grid)
mx, my = size(config.cells)
xx, yy = config.grid.x, config.grid.y
δ = ones(nx, ny)
εr  = ones(mx, my, 1)
bcs = zeros(Int8, nx, ny, 1)
bcs[ nx,  1, 1] = 1
bcs[ nx, ny, 1] = 2
set_permittivity(εr)
add_electrode(bcs .== 1, +1e3V)
add_electrode(bcs .== 2, -1e3V)
############################################
import ParticleInCell
import Diagnostics

function ParticleInCell.enter_loop()
  Diagnostics.open_container("06-field")
  Diagnostics.open_container("06-particle")
end

function ParticleInCell.after_loop(it)
  Diagnostics.save_diagnostic("E",   "06-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("ϕ",   "06-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("ν",   "06-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("nO",  "06-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("ne-", "06-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("nO+", "06-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("pve-","06-particle",it, Δt*it-Δt)
  Diagnostics.save_diagnostic("pvO+","06-particle",it, Δt*it-Δt)
end

function ParticleInCell.exit_loop()
  Diagnostics.close_container("06-field")
  Diagnostics.close_container("06-particle")
end
ParticleInCell.init(ParticleInCell.MaxwellianSource(5e3/Δt, [0 Lx; 0 Ly], [.5e6 -1e6; .5e6 -1e6]), e, Δt)
ParticleInCell.init(ParticleInCell.DensitySource(5e3δ, config.grid), O, Δt)
@time ParticleInCell.solve(config, Δt, ts, ε0)