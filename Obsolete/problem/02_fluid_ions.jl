# + spatial and temporal paramters
nx = 20         # number of nodes in x direction
ny = 20         # number of nodes in y direction
ts = 1000       # number of time steps
Δh = 5cm        # cell size
Δt = 5000ns     # time step
Lx = nx*Δh      # domain length in x direction
Ly = ny*Δh      # domain length in y direction

# + species and sources
xs = 0m:Δh:Lx
ys = 0m:Δh:Ly
iO = create_fluid_species("O+", 1.0, +1qe, 8mp, nx+1, ny+1)

# + grid, solver and pusher
using RegularGrids, FiniteDifferenceMethod, ParticleInCell
grid    = create_uniform_grid(xs, ys)
pusher  = nothing
solver  = create_poisson_solver(grid, ε0)
species = [iO]
interactions = []

# + boundary conditions
nx, ny = size(grid)
xx, yy = grid.coords
δ = @. exp(-(xx-0.5Lx)^2/0.03Lx -(yy-0.5Ly)^2/0.03Ly)
bcs = zeros(Int8, nx, ny)
bcs[ nx,  1] = 1
bcs[ nx, ny] = 2
create_electrode(bcs .== 1, solver, grid; σ=30ε0)
create_electrode(bcs .== 2, solver, grid; fixed=true)

function ParticleInCell.after_push(part, grid)
 ParticleInCell.wrap!(part, grid)
end

# + hooks
function start(dt)
  iO.n .= 0.0
  init(DensitySource(1e5δ, grid), iO, Δt)
end

using Diagnostics
function iteration(i, t, dt)
  cd("/tmp")
  new_iteration("02_fluid_ions", i, t, dt) do it
    save_record(it, "nO+")
    save_record(it, "dO+")
    save_record(it, "vO+")
    save_record(it, "phi")
    save_record(it, "E")
  end
  [(:iteration, i)]
end

using XDMF
function stop()
  println("Exporting to XDMF...")
  cd("/tmp/02_fluid_ions")
  fields = new_document()
  xdmf(1:ts) do it
    write_fields(it, fields)
  end
  save_document(fields, "fields")
end