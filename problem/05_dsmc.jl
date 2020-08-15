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
xx, yy = config.grid.coords
bcs = zeros(Int8, nx, ny)
bcs[ nx,  1] = 1
bcs[ nx, ny] = 2
create_electrode(bcs .== 1, config; σ=1e3ε0)
create_electrode(bcs .== 2, config; fixed=true)
############################################
import ParticleInCell
using Diagnostics
using XDMF

function thermal_speed(T, m)
  kB = 1.3806503e-23
  sqrt(2kB*T/m)
end

function ParticleInCell.after_loop(i, t, dt)
  cd("/tmp")
  new_iteration("05_dsmc", i, t, dt) do it
    save_record(it, "phi")
    save_record(it, "nuDSMC")
    save_record(it, "nO")
    save_record(it, "ne-")
    save_record(it, "nO+")
    save_record(it, "E")
    save_records(it, "e-/")
    save_records(it, "O+/")
    save_records(it, "O/")
  end
end

function ParticleInCell.exit_loop()
  println("Exporting to XDMF...")
  cd("/tmp/05_dsmc")
  electrons = new_document()
  neutrals = new_document()
  fields = new_document()
  ions = new_document()
  xdmf(1:ts) do it
    write_species(it, electrons, "e-")
    write_species(it, neutrals, "O")
    write_species(it, ions, "O+")
    write_fields(it, fields)
  end
  save_document(electrons, "electrons")
  save_document(neutrals, "neutrals")
  save_document(fields, "fields")
  save_document(ions, "ions")
end

νth = thermal_speed(300, O.m)
ParticleInCell.init(ParticleInCell.MaxwellianSource{2,3}(5e3/Δt, [1.0Lx 1.0Ly], [.5e6 .5e6 0]), e, Δt)
ParticleInCell.init(ParticleInCell.MaxwellianSource{2,3}(5e3/Δt, [1.0Lx 1.0Ly], [νth νth νth]), O, Δt)
@time ParticleInCell.solve(config, Δt, ts)
