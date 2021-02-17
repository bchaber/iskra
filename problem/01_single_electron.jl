# + spatial and temporal paramters
nx = 20         # number of nodes in x direction
ny = 20         # number of nodes in y direction
ts = 1000       # number of time steps
Δh = 5cm        # cell size
Δt = 20ns       # time step
Lx = nx*Δh      # domain length in x direction
Ly = ny*Δh      # domain length in y direction

# + species and sources
xs = 0m:Δh:Lx
ys = 0m:Δh:Ly
e = create_kinetic_species("e-", 20_000,-1qe, 1me, 1);
γ = create_thermalized_beam(e, [Lx Ly], [+0.05Δh/Δt 0. 0.]; T=300K, rate=1.0/Δt)

# + grid, solver and pusher
using RegularGrids, FiniteDifferenceMethod, ParticleInCell
grid    = create_uniform_grid(xs, ys)
solver  = create_poisson_solver(grid, ε0)
pusher  = create_boris_pusher()
species = [e]
interactions = []

# + boundary conditions
nx, ny = size(grid)
bcs = zeros(Int8, nx, ny)
bcs[ nx, ny] = 1
bcs[ nx,  1] = 2

create_electrode(bcs .== 1, solver, grid; σ=-1ε0)
create_electrode(bcs .== 2, solver, grid; fixed=true)

function ParticleInCell.after_push(part, grid)
 ParticleInCell.discard!(part, grid)
end

# + hooks
function start(dt)
  e.np = 0
  init(γ, e, dt)
end

using Diagnostics
function iteration(i, t, dt)
  cd("/tmp")
  new_iteration("01_single_particle", i, t, dt) do it
    save_records(it, "e-/")
    save_record(it, "rho")
    save_record(it, "phi")
  end
  [(:iteration, i), (:e, e.np)]
end

using XDMF
function stop()
  println("Exporting to XDMF...")
  cd("/tmp/01_single_particle")
  electrons = new_document()
  fields = new_document()
  xdmf(1:ts) do it
    write_species(it, electrons, "e-")
    write_fields(it, fields)
  end
  save_document(electrons, "electrons")
  save_document(fields, "fields")
end
#+
