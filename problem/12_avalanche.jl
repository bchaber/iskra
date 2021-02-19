# + spatial and temporal paramters
Td  = 1e-21u"V/m"
nAr = 1e22 #41.4u"Pa" / (1.380649e-23u"J/K" * 300u"K")
T  = 300.0K
E  = 5000 #500Td * nAr
d  = 0.04 #< 500u"V"/E
nx = 32
ny = 64
ts = 1024
Δh = d/nx

simulationVolume = nx * Δh * ny * Δh
numCells         = nx * ny

Δt = 0.075ns
Lx = nx*Δh
Ly = ny*Δh
println("nAr: ", nAr, " d: ", d)
println("Δt: ", Δt, "\nΔh: ", Δh, "\nLx: ", Lx)

# + species and sources
xs = 0m:Δh:Lx
ys = 0m:Δh:Ly
e  = create_kinetic_species("e-",   100_000,-1qe, 1.00me, 1);
iAr = create_kinetic_species("Ar+", 100_000,+1qe, 3.99mp, 1);
Ar = FluidSpecies("Ar", 1.0, 0qe, 3.99mp, nAr*ones(nx+1, ny+1), T)
se   = create_thermalized_beam(e,   [Δh Δh], [0. 0. 0.]; dx=[0.0 Ly/2], T=T, rate=1/Δt)

# + grid, solver and pusher
using Chemistry
include("../Chemistry/src/biagi71-Ar.jl")
electron = mcc(@reactions begin
    σ₁, e + Ar --> e + Ar
    σ₂, e + Ar --> e + Ar, Chemistry.MCC.Excitation(11.55eV)
    σ₃, e + Ar --> e + Ar, Chemistry.MCC.Excitation(13.00eV)
    σ₄, e + Ar --> e + e + iAr, Chemistry.MCC.Ionization(15.7eV)
end)

using RegularGrids, FiniteDifferenceMethod, ParticleInCell
grid    = create_uniform_grid(xs, ys)
solver  = create_poisson_solver(grid, ε0)
pusher  = create_boris_pusher()
species = [e, iAr, Ar]
interactions = [electron]

# + boundary conditions
nx, ny = size(grid)
bcs = zeros(Int8, nx, ny)
bcs[ 1, :] .= 1
bcs[nx, :] .= 2
apply_periodic(solver, 1)
apply_dirichlet(solver, bcs .== 1, 0.0)
apply_dirichlet(solver, bcs .== 2, E*d)

function ParticleInCell.after_push(part, grid)
  discarded = 
    ParticleInCell.discard!(part, grid; dims=1)
  ParticleInCell.wrap!(part, grid; dims=2)
end

# + hooks
function start(dt)
  init(se, e, Δt)
end

using Diagnostics
function iteration(i, t, dt)
  cd("/tmp")
  new_iteration("12_avalanche", i, t, dt) do it
    save_records(it, "e-/")
    save_records(it, "Ar+/")
    save_record(it, "nuMCC-e--1")
    save_record(it, "nuMCC-e--2")
    save_record(it, "nuMCC-e--3")
    save_record(it, "nuMCC-e--4")
    save_record(it, "rho")
    save_record(it, "phi")
    save_record(it, "E")
    save_record(it, "nAr+")
    save_record(it, "ne-")
  end
  [(:iteration, i), (:electrons, e.np), (:ions, iAr.np)]
end

using XDMF
function stop()
  println("Exporting to XDMF...")
  cd("/tmp/12_avalanche")
  electrons = new_document()
  ions = new_document()
  fields = new_document()
  xdmf(1:ts) do it
    write_species(it, electrons, "e-")
    write_species(it, ions, "Ar+")
    write_fields(it, fields)
  end
  save_document(electrons, "electrons")
  save_document(ions, "ions")
  save_document(fields, "fields")
end