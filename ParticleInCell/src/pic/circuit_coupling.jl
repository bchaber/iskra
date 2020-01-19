using Circuit
using FiniteDifferenceMethod

const DegreeOfFreedom = SubArray{Float64,0,Array{Float64,1},Tuple{Int64},true}
struct FixedPotentialElectrode <: Surface
  ϕ :: DegreeOfFreedom
  dq :: Float64
  area :: Float64
end

mutable struct FloatingPotentialElectrode <: Surface
  σ :: DegreeOfFreedom
  ϕ :: DegreeOfFreedom
  dq :: Float64
  area :: Float64
end

struct PlasmaDevice
  fixed :: Vector{FixedPotentialElectrode}
  floating :: Vector{FloatingPotentialElectrode}
end

struct PlasmaCircuit
  pd :: PlasmaDevice
end

PlasmaDevice() =
  PlasmaDevice(FixedPotentialElectrode[], FloatingPotentialElectrode[])
function assign!(pd :: PlasmaDevice, s :: Surface) end
function assign!(pd :: PlasmaDevice, s :: FixedPotentialElectrode)
  push!(pd.fixed, s)
end
function assign!(pd :: PlasmaDevice, s :: FloatingPotentialElectrode)
  push!(pd.floating, s)
end

function create_plasma_circuit(tracker :: SurfaceTracker)
  pd = PlasmaDevice()
  for surface in tracker.surfaces
  	assign!(pd, surface)
  end

  if length(pd.floating) > 0 && length(pd.fixed) < 1
  	println("At least one electrode should have fixed potential!")
  end

  PlasmaCircuit(pd)
end

function advance!(circuit :: Nothing, ϕ, Δt, config, ε0) end
function advance!(circuit :: PlasmaCircuit, ϕ, Δt, config, ε0)
  for electrode in circuit.pd.floating
  	electrode.σ .+= electrode.dq/ε0 # TODO
  	electrode.dq = 0.0
  end
end

function advance!(circuit :: CircuitRLC, ϕ, Δt, config, ε0)
  nx, ny = size(config.grid)
  Δx, Δy, Δz = config.grid.Δh
  A = ny * Δy*Δz
  V = ϕ[1,1] - ϕ[nx,ny]
  advance_circuit!(circuit, V, Δt)
  dσ = -Δt*circuit.i/A
  @diag "dσ" TimeData(dσ)
  σ = FiniteDifferenceMethod.get_rhs(config.solver, :σ, 1)
  σ .+= dσ/ε0
end

function hit!(s::FloatingPotentialElectrode,
              part::KineticSpecies,
              st ::SurfaceTracker,
              pt ::TrackedParticle,
              pt′::TrackedParticle)
  _, p  = pt
  dq    = part.q*part.wg[p]
  s.dq += dq
  s.σ .+= dq/s.area
  absorbed!(st, pt)
end

function hit!(s::FixedPotentialElectrode,
              part::KineticSpecies,
              st ::SurfaceTracker,
              pt ::TrackedParticle,
              pt′::TrackedParticle)
  absorbed!(st, pt)
end