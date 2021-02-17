# + spatial and temporal paramters
nx = 20         # number of nodes in x direction
ny = 20         # number of nodes in y direction
ts = 1000       # number of time steps
Δh = 5cm        # cell size
Δt = 50ns       # time step
Lx = nx*Δh      # domain length in x direction
Ly = ny*Δh      # domain length in y direction

# + species and sources
xs = 0m:Δh:Lx
ys = 0m:Δh:Ly
O  = create_fluid_species("O",  0.0, 0qe, 8mp, nx+1, ny+1)
iO = create_fluid_species("O+", 0.0,+1qe, 8mp, nx+1, ny+1)
e  = create_fluid_species("e-", 1.0,-1qe, 1me, nx+1, ny+1)

# + grid, solver and pusher
using Chemistry
using RegularGrid, FiniteDifferenceMethod, ParticleInCell
grid    = create_uniform_grid(xs, ys)
solver  = create_poisson_solver(grid, ε0)
pusher  = create_boris_pusher()
species = [O, iO, e]

#σ = CrossSection(0:0.3:1.5, [0, 0.1e-4, 0.4e-4, 0.5e-4, 0.7e-4, 0.9e-4])
σ(E) = 2
chemistry = chemical(@reactions begin
    σ, e + O --> 2e + iO
end)
interactions = [chemistry]

# + boundary conditions
nx, ny = size(grid)
xx, yy = grid.coords
δ = @. exp(-(xx-0.5Lx)^2/0.03Lx -(yy-0.5Ly)^2/0.03Ly)
bcs = zeros(Int8, nx, ny)
bcs[ nx,  1] = 1
bcs[ nx, ny] = 2
create_electrode(bcs .== 1, solver, grid; σ=1ε0)
create_electrode(bcs .== 2, solver, grid; fixed=true)

function ParticleInCell.after_push(part, grid)
 ParticleInCell.wrap!(part, grid)
end

# + hooks
function start(dt)
  O.n .= 0.0
  e.n .= 0.0  
  init(DensitySource(1e6δ, grid), O, Δt)
  init(DensitySource(1e4δ, grid), e, Δt)
end

using Diagnostics
function iteration(i, t, dt)
  cd("/tmp")
  new_iteration("03_chemistry", i, t, dt) do it
    save_record(it, "phi")
    save_record(it, "nO")
    save_record(it, "nO+")
    save_record(it, "ne-")
    save_record(it, "E")
  end
  [(:iteration, i)]
end

using XDMF
function stop()
  println("Exporting to XDMF...")
  cd("/tmp/03_chemistry")
  fields = new_document()
  xdmf(1:ts) do it
    write_fields(it, fields)
  end
  save_document(fields, "fields")
end