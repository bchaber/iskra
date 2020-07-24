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

struct PlasmaDevice <: Circuit.CircuitDevice
  positive :: FloatingPotentialElectrode
  negative :: FixedPotentialElectrode
end

function Circuit.voltage(pd :: PlasmaDevice)
  pd.positive.ϕ - pd.negative.ϕ
end
function foo!(pd :: Circuit.CircuitDevice, current, Δt)
  0.0
end
function foo!(pd :: PlasmaDevice, current, Δt)
  positive, negative = pd.positive, pd.negative
  positive.dq = .0
  dσ = -Δt*current/positive.area
end
function advance!(circuit :: Nothing, ϕ, Δt, config) end
function advance!(circuit :: CircuitRLC, ϕ, Δt, config)
  advance_circuit!(circuit, 0, Δt)
  dσ = foo!(circuit.ext, circuit.i, Δt)
  ε0 = config.solver.ε0
  σ = get_rhs(config.solver, :σ, 1) # Hack!
  σ .+= dσ
  @probe "dsigma" dσ "C/m^2"
  @probe  "sigma" first(σ) "C/m^2"
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