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
Δt = 1ns        # time step
Lx = nx*Δh      # domain length in x direction
Ly = ny*Δh      # domain length in y direction
############################################
xs, ys = 0m:Δh:Lx, 0m:Δh:Ly
sx, sv = [0 Lx; 0 Ly], [0 0; 0 0]
O  = create_kinetic_species("O",  20_000, 0qe, 8mp, 1)
e  = create_kinetic_species("e-", 20_000,-1qe, 1me, 1)
iO = create_kinetic_species("O+", 20_000,+1qe, 8mp, 1)
using Chemistry
import RegularGrid, FiniteDifferenceMethod, ParticleInCell
config.grid    = RegularGrid.create_uniform_grid(xs, ys)
config.cells   = RegularGrid.create_staggered_grid(config.grid)
config.solver  = FiniteDifferenceMethod.create_poisson_solver(config.grid, ε0)
config.pusher  = ParticleInCell.create_boris_pusher()
config.species = [e, O, iO]

σ = CrossSection(3e6:1e6:6e6, [0.01, 0.1, 2.0, 0.01])
collisions = dsmc(@reactions begin
    σ, e + O --> O + e
end)
config.interactions = [collisions]
############################################
nx, ny = size(config.grid)
mx, my = size(config.cells)
xx, yy = config.grid.x, config.grid.y
εr  = ones(mx, my, 1)
bcs = zeros(Int8, nx, ny, 1)
bcs[ nx,  1, 1] = 1
bcs[ nx, ny, 1] = 2
set_permittivity(εr)
create_electrode(bcs .== 1, config; σ=1e3ε0)
create_electrode(bcs .== 2, config; fixed=true)
############################################
import ParticleInCell
import Diagnostics

function thermal_speed(T, m)
  kB = 1.3806503e-23
  sqrt(2kB*T/m)
end

function ParticleInCell.enter_loop()
  Diagnostics.open_container("05-field")
  Diagnostics.open_container("05-particle")
end

function ParticleInCell.after_loop(it)
  Diagnostics.save_diagnostic("ν",   "05-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("E",   "05-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("ϕ",   "05-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("nO",  "05-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("ne-", "05-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("nO+", "05-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("pvO", "05-particle",it, Δt*it-Δt)
  Diagnostics.save_diagnostic("pve-","05-particle",it, Δt*it-Δt)
  Diagnostics.save_diagnostic("pvO+","05-particle",it, Δt*it-Δt)
end

function ParticleInCell.exit_loop()
  Diagnostics.close_container("05-field")
  Diagnostics.close_container("05-particle")
end
νth = thermal_speed(300, O.m)
ParticleInCell.init(ParticleInCell.MaxwellianSource(5e3/Δt, [0 Lx; 0 Ly], [.5e6 -1e6; .5e6 -1e6]), e, Δt)
ParticleInCell.init(ParticleInCell.MaxwellianSource(5e3/Δt, [0 Lx; 0 Ly], [0.     νth; 0.   νth]), O, Δt)
@time ParticleInCell.solve(config, Δt, ts)
