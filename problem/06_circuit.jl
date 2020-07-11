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
Δt = 10ns       # time step
Lx = nx*Δh      # domain length in x direction
Ly = ny*Δh      # domain length in y direction
############################################
xs, ys = 0m:Δh:Lx, 0m:Δh:Ly
sx, sv = [0 Lx; 0 Ly], [0 0; 0 0]
O  = create_fluid_species("O", 1.0, 0qe, 8mp, nx+1, ny+1)
e  = create_kinetic_species("e-", 20_000,-1qe, 1me, 50e3)
using Chemistry, Circuit
import RegularGrid, FiniteDifferenceMethod, ParticleInCell
config.grid    = RegularGrid.create_uniform_grid(xs, ys)
config.cells   = RegularGrid.create_staggered_grid(config.grid)
config.solver  = FiniteDifferenceMethod.create_poisson_solver(config.grid, ε0)
config.pusher  = ParticleInCell.create_boris_pusher()
config.species = [e, O]
############################################
nx, ny = size(config.grid)
mx, my = size(config.cells)
xx, yy = config.grid.coords
δ = ones(nx, ny)
εr  = ones(mx, my, 1)
bcs = zeros(Int8, nx, ny, 1)
bcs[ 1, 1:ny, 1] .= 1
bcs[nx, 1:ny, 1] .= 2
set_permittivity(εr)
driven   = create_electrode(bcs .== 1, config; σ=1ε0)
grounded = create_electrode(bcs .== 2, config; fixed=true)

config.circuit = rlc(@netlist begin
  V1, 3, GND, t -> sin(2π*.1e3*t)
  L1, NOD, VCC, 1nH
  R1, GND, NOD, 1
end)
############################################
import ParticleInCell
using Diagnostics
using XDMF

function ParticleInCell.after_loop(i, t, dt)
  cd("/tmp")
  new_iteration("06_circuit", i, t, dt) do it
    save_records(it, "e-/")
    #save_record(it, "rho")
    #save_record(it, "phi")
    #save_record(it, "ne-")
    #save_record(it, "nO")
    #save_record(it, "E")
    save_record(it, "Q1")
    save_record(it, "I1")
    save_record(it, "V1")
    save_record(it, "Vext")
  end
end

function ParticleInCell.exit_loop()
  println("Exporting to XDMF...")
  cd("/tmp/06_circuit")
  electrons = new_document()
  probes = new_document()
  fields = new_document()
  xdmf(1:ts) do it
    write_species(it, electrons, "e-")
    write_probes(it, probes)
  end
  save_document(electrons, "electrons")
  save_document(probes, "probes")
end

ParticleInCell.init(ParticleInCell.MaxwellianSource(1e3/Δt, [0 Lx; 0 Ly], [0 0; 0 0]), e, Δt)
ParticleInCell.init(ParticleInCell.DensitySource(0δ, config.grid), O, Δt)
@time ParticleInCell.solve(config, Δt, ts)