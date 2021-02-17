# + spatial and temporal paramters
nx = 20         # number of nodes in x direction
ny = 20         # number of nodes in y direction
ts = 1000       # number of time steps
Δh = 5cm        # cell size
Δt = 1ns        # time step
Lx = nx*Δh      # domain length in x direction
Ly = ny*Δh      # domain length in y direction

# + species and sources
xs = 0m:Δh:Lx
ys = 0m:Δh:Ly
O  = create_kinetic_species("O",  20_000, 0qe, 8mp, 1)
e  = create_kinetic_species("e-", 20_000,-1qe, 1me, 1)
iO = create_kinetic_species("O+", 20_000,+1qe, 8mp, 1)

# + grid, solver and pusher
using Chemistry
using RegularGrid, FiniteDifferenceMethod, ParticleInCell
grid    = create_uniform_grid(xs, ys)
solver  = create_poisson_solver(grid, ε0)
pusher  = create_boris_pusher()
species = [e, O, iO]

σ = CrossSection(3e6:1e6:6e6, [0.01, 0.1, 2.0, 0.01])
collisions = dsmc(@reactions begin
    σ, e + O --> O + e
end)
interactions = [collisions]

# + boundary conditions
nx, ny = size(grid)
xx, yy = grid.coords
bcs = zeros(Int8, nx, ny)
bcs[ nx,  1] = 1
bcs[ nx, ny] = 2
create_electrode(bcs .== 1, solver, grid; σ=1e3ε0)
create_electrode(bcs .== 2, solver, grid; fixed=true)

function ParticleInCell.after_push(part, grid)
 ParticleInCell.wrap!(part, grid)
end

# + hooks
function thermal_speed(T, m)
  kB = 1.3806503e-23
  sqrt(2kB*T/m)
end

function start(dt)
  e.np = 0
  O.np = 0
  νth = thermal_speed(300, O.m)
  init(MaxwellianSource{2,3}(5e3/Δt, [1.0Lx 1.0Ly], [.5e6 .5e6 0]), e, Δt)
  init(MaxwellianSource{2,3}(5e3/Δt, [1.0Lx 1.0Ly], [νth νth νth]), O, Δt)
end

using Diagnostics
function iteration(i, t, dt)
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
  [(:iteration, i)]
end

using XDMF
function stop()
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