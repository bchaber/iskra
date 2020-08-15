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
Δt = 50ns       # time step
Lx = nx*Δh      # domain length in x direction
Ly = ny*Δh      # domain length in y direction
############################################
xs, ys = 0m:Δh:Lx, 0m:Δh:Ly
O  = create_fluid_species("O",  0.0, 0qe, 8mp, nx+1, ny+1)
iO = create_fluid_species("O+", 0.0,+1qe, 8mp, nx+1, ny+1)
e  = create_fluid_species("e-", 1.0,-1qe, 1me, nx+1, ny+1)
using Chemistry
import RegularGrid, FiniteDifferenceMethod, ParticleInCell
config.grid    = RegularGrid.create_uniform_grid(xs, ys)
config.cells   = RegularGrid.create_staggered_grid(config.grid)
config.solver  = FiniteDifferenceMethod.create_poisson_solver(config.grid, ε0)
config.pusher  = ParticleInCell.create_boris_pusher()
config.species = [O, iO, e]

#σ = CrossSection(0:0.3:1.5, [0, 0.1e-4, 0.4e-4, 0.5e-4, 0.7e-4, 0.9e-4])
σ(E) = 2
chemistry = chemical(@reactions begin
    σ, e + O --> 2e + iO
end)
config.interactions = [chemistry]

############################################
nx, ny = size(config.grid)
mx, my = size(config.cells)
xx, yy = config.grid.coords
δ = @. exp(-(xx-0.5Lx)^2/0.03Lx -(yy-0.5Ly)^2/0.03Ly)
bcs = zeros(Int8, nx, ny)
bcs[ nx,  1] = 1
bcs[ nx, ny] = 2
create_electrode(bcs .== 1, config; σ=1ε0)
create_electrode(bcs .== 2, config; fixed=true)
############################################
import ParticleInCell
using Diagnostics
using XDMF

function ParticleInCell.after_loop(i, t, dt)
  cd("/tmp")
  new_iteration("03_chemistry", i, t, dt) do it
    save_record(it, "phi")
    save_record(it, "nO")
    save_record(it, "nO+")
    save_record(it, "ne-")
    save_record(it, "E")
  end
end

function ParticleInCell.exit_loop()
  println("Exporting to XDMF...")
  cd("/tmp/03_chemistry")
  fields = new_document()
  xdmf(1:ts) do it
    write_fields(it, fields)
  end
  save_document(fields, "fields")
end

ParticleInCell.init(ParticleInCell.DensitySource(1e6δ, config.grid), O, Δt)
ParticleInCell.init(ParticleInCell.DensitySource(1e4δ, config.grid), e, Δt)
@time ParticleInCell.solve(config, Δt, ts)