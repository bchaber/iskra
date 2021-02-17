# Turner Benchmark
#

# + spatial and temporal paramters
Te = 30_000K
Ti = 300K
Tn = 300K
L  = 6.7cm
nHe = 9.64e20
ne = 2.56e14
f  = 13.56MHz
nx = 128
ny = 1
ds = L/nx
ts = 1_000
Δh = ds

simulationVolume = nx * Δh * ny * Δh
numCells         = nx * ny

electronParticles    = 128 * numCells
totalNumElectrons    = ne * simulationVolume
electronNumRatio     = totalNumElectrons / electronParticles

ionParticles         = 128 * numCells
totalNumIons         = nHe * simulationVolume
ionNumRatio          = totalNumIons / ionParticles

Δt = 1/400f
Lx = nx*Δh
Ly = ny*Δh

println("Δt: ", Δt, "\nΔh: ", Δh, "\nLx: ", Lx)
println("electrons: ", electronParticles, " wg: ", electronNumRatio)
println("ions: ", ionParticles, " wg: ", ionNumRatio)
println("kB Te / me: ", .5thermal_speed(Te, me))

# + species and sources
xs = 0m:Δh:Lx
ys = 0m:Δh:Ly
e  = create_kinetic_species("e-", 50_000,-1qe, 1me, electronNumRatio);
iHe = create_kinetic_species("He+", 50_000,+1qe, 3.99mp, electronNumRatio);
He  = FluidSpecies("He", 1.0, 0qe, 3.99mp, nHe*ones(nx+1, ny+1), float(Tn))

se  = create_thermalized_beam(e,  [Lx Ly], [0. 0. 0.]; T=Te, rate=electronParticles/Δt)
siHe = create_thermalized_beam(iHe, [Lx Ly], [0. 0. 0.]; T=Ti, rate=ionParticles/Δt)

# + grid, solver and pusher
using Chemistry
include("../Chemistry/src/biagi71-He.jl")
include("../Chemistry/src/phelps-He.jl")
electron = mcc(@reactions begin
    σ₁, e + He --> e + He
    σ₂, e + He --> e + He, Chemistry.MCC.Excitation(19.82eV)
    σ₃, e + He --> e + He, Chemistry.MCC.Excitation(20.61eV)
    σ₄, e + He --> e + e + iHe, Chemistry.MCC.Ionization(24.587eV)
end)

ion = mcc(@reactions begin
    σₑ₂, iHe + He --> iHe + He, Chemistry.MCC.ElasticBackward()
    σₑ₁, iHe + He --> iHe + He, Chemistry.MCC.ElasticIsotropic()
end)

using RegularGrids, FiniteDifferenceMethod, ParticleInCell
grid    = create_uniform_grid(xs, ys)
solver  = create_poisson_solver(grid, ε0)
pusher  = create_boris_pusher()
species = [e, iHe, He]
interactions = [electron, ion]

# + boundary conditions
nx, ny = size(grid)
bcs = zeros(Int8, nx, ny)
bcs[ 1, :] .= 1
bcs[nx, :] .= 2
apply_periodic(solver, 1)
apply_dirichlet(solver, bcs .== 1, 0.0)
apply_dirichlet(solver, bcs .== 2, 0.0)

function ParticleInCell.after_push(part, grid)
  ParticleInCell.discard!(part, grid; dims=1)
  ParticleInCell.wrap!(part, grid; dims=2)
end

# + hooks
function start(dt)
  e.np = 0
  iHe.np = 0
  init(se,  e, Δt)
  init(siHe,iHe, Δt)
end

using Diagnostics
function iteration(i, t, dt)
  apply_dirichlet(solver, bcs .== 1, 450sin(2π*f*t))
  if i % 1 != 1
    cd("/tmp")
    new_iteration("11_rf_discharge", i, t, dt) do it
      save_records(it, "e-/")
      save_records(it, "He+/")
      save_record(it, "nuMCC-e--1")
      save_record(it, "nuMCC-e--2")
      save_record(it, "nuMCC-e--3")
      save_record(it, "nuMCC-e--4")
      save_record(it, "nuMCC-He+-1")
      save_record(it, "nuMCC-He+-2")
      save_record(it, "rho")
      save_record(it, "phi")
      save_record(it, "E")
      save_record(it, "nHe+")
      save_record(it, "ne-")
    end
  end
  [(:iteration, i), (:e, e.np), (:iHe, iHe.np)]
end

using XDMF
function stop()
end