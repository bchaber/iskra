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
ts = 250        # number of time steps
Δh = 5cm        # cell size
Δt = 0.001ns    # time step
Lx = nx*Δh      # domain length in x direction
Ly = ny*Δh      # domain length in y direction
############################################
xs, ys = 0m:Δh:Lx, 0m:Δh:Ly
sx, sv = [0 Lx; 0 Ly], [0 0; 0 0]
O  = create_fluid_species("O", 1.0, 0qe, 8mp, nx+1, ny+1)
e  = create_kinetic_species("e-", 20_000,-1qe, 1me, 1)
iO = create_kinetic_species("O+", 20_000,+1qe, 8mp/2000, 1)
using Chemistry, Circuit
import RegularGrid, FiniteDifferenceMethod, ParticleInCell
config.grid    = RegularGrid.create_uniform_grid(xs, ys)
config.cells   = RegularGrid.create_staggered_grid(config.grid)
config.solver  = FiniteDifferenceMethod.create_poisson_solver(config.grid)
config.pusher  = ParticleInCell.create_boris_pusher()
config.species = [e, O, iO]
config.circuit = rlc(@netlist begin
  V1, 3, GND, t -> sin(2π*.1e6*t)
  L1, NOD, VCC, 1e-12
  R1, GND, NOD, 1e-6
end)
############################################
nx, ny = size(config.grid)
mx, my = size(config.cells)
xx, yy = config.grid.x, config.grid.y
δ = ones(nx, ny)
εr  = ones(mx, my, 1)
bcs = zeros(Int8, nx, ny, 1)
bcs[ 1, 1:ny, 1] .= 1
bcs[nx, 1:ny, 1] .= 2
set_permittivity(εr)
FiniteDifferenceMethod.add_new_dof(config.solver, :σ)
FiniteDifferenceMethod.get_rhs(config.solver, :σ, 1) .= 1.
FiniteDifferenceMethod.apply_neumann(config.solver, bcs .== 1, 1)
FiniteDifferenceMethod.apply_dirichlet(config.solver, bcs .== 2, 0)
############################################
import ParticleInCell
import Diagnostics

function ParticleInCell.enter_loop()
  Diagnostics.open_container("06-field")
  Diagnostics.open_container("06-particle")
  Diagnostics.open_container("06-circuit")
end

function ParticleInCell.after_loop(it)
  Diagnostics.save_diagnostic("E",   "06-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("ϕ",   "06-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("nO",  "06-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("ne-", "06-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("nO+", "06-field",   it, Δt*it-Δt)
  Diagnostics.save_diagnostic("pve-","06-particle",it, Δt*it-Δt)
  Diagnostics.save_diagnostic("pvO+","06-particle",it, Δt*it-Δt)
  Diagnostics.save_diagnostic("i",   "06-circuit", it, Δt*it-Δt)
  Diagnostics.save_diagnostic("q",   "06-circuit", it, Δt*it-Δt)
  Diagnostics.save_diagnostic("V",   "06-circuit", it, Δt*it-Δt)
  Diagnostics.save_diagnostic("dσ",  "06-circuit", it, Δt*it-Δt)
  Diagnostics.save_diagnostic("Vext","06-circuit", it, Δt*it-Δt)
end

function ParticleInCell.exit_loop()
  Diagnostics.close_container("06-field")
  Diagnostics.close_container("06-particle")
  Diagnostics.close_container("06-circuit")
end
ParticleInCell.init(ParticleInCell.MaxwellianSource(1e3/Δt, [0 Lx; 0 Ly], [0 0; 0 0]), iO, Δt)
ParticleInCell.init(ParticleInCell.DensitySource(0δ, config.grid), O, Δt)
@time ParticleInCell.solve(config, Δt, ts, ε0)