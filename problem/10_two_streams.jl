# Two Stream Instability
# Simulated in quasi-2D with only one cell in y-axis.
#
# Parameters are based on the simulation of two_stream_ee_es.inp distributed with XOOPIC.
# Two electron beams with opposite velocities interact with each other. There is also a fluid species of ions that neutralizes electrons' charge. In XOOPIC simulation ions (He+) are also a kinetic species.

# + spatial and temporal paramters
νtherm = thermal_speed(300K, me)
νdrift = 1e7 #m/s
mHe = 4.002602me/5.48579903e-04
nHe = 1e24
ne = 1nHe
f  = 9kHz
ω  = 2π * f * sqrt(2e-6ne)
ds = 5e-3c0/ω
ts = 1000
nx = 128
ny = 1
Δh = ds

simulationVolume = nx * Δh * ny * Δh
numCells         = nx * ny
electronParticles    = 10 * numCells
totalNumElectrons    = ne * simulationVolume
electronNumRatio     = totalNumElectrons / electronParticles
Δt = .4Δh/νdrift/sqrt(2.0)
Lx = nx*Δh
Ly = ny*Δh

println("Δt: ", Δt, "\nΔh: ", Δh, "\nLx: ", Lx)
println("νtherm: ", round(νtherm/c0; sigdigits=1), "c")
println("νdrift: ", round(νdrift/c0; sigdigits=3), "c")
println("electrons: ", electronParticles, " wg: ", electronNumRatio)

# + species and sources
xs = 0m:Δh:Lx
ys = 0m:Δh:Ly
e   = create_kinetic_species("e-", 20_000,-1qe, 1me, electronNumRatio);
iHe = create_kinetic_species("He+", 20_000, +1qe, 1mHe, electronNumRatio)

fwd = create_thermalized_beam(e, [Lx Ly], [+νdrift 0 0]; T=300.0K, rate=electronParticles/2/Δt)
rev = create_thermalized_beam(e, [Lx Ly], [-νdrift 0 0]; T=300.0K, rate=electronParticles/2/Δt)

# + grid, solver and pusher
using RegularGrid, FiniteDifferenceMethod, ParticleInCell
grid    = create_uniform_grid(xs, ys)
solver  = create_poisson_solver(grid, ε0)
pusher  = create_boris_pusher()
species = [e, iHe]
interactions = []

# + boundary conditions
apply_periodic(solver, 1)
apply_periodic(solver, 2)

function ParticleInCell.after_push(part, grid)
 ParticleInCell.wrap!(part, grid)
end

# + hooks
function start(dt)
  init(fwd, e, Δt)
  init(rev, e, Δt)
  # place ions in the same positions as electrons
  iHe.x .= e.x
  iHe.v .= 
  iHe.np = e.np
end

using Diagnostics
function iteration(i, t, dt)
  cd("/tmp")
  new_iteration("10_two_streams", i, t, dt) do it
    save_records(it, "e-/")
    save_records(it, "He+/")
    save_record(it, "rho")
    save_record(it, "phi")
    save_record(it, "E")
    save_record(it, "nHe+")
    save_record(it, "ne-")
  end
  [(:iteration, i), (:e, e.np)]
end

using XDMF
function stop()
  println("Exporting to XDMF...")
  cd("/tmp/10_two_streams")
  electrons = new_document()
  fields = new_document()
  xdmf(1:ts) do it
    write_species(it, electrons, "e-")
    write_fields(it, fields)
  end
  save_document(electrons, "electrons")
  save_document(fields, "fields")
end
