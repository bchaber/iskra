module Circuit
using DataStructures

export CircuitRLC
export @netlist
export rlc
export advance_circuit!
export TimeData
export resonant_frequency, damping_factor

import Diagnostics: DiagnosticData, @diag, save_diagnostic
import PlotVTK: pvd_add_timestep, field_as_vectors
struct TimeData <: DiagnosticData
   y :: Array{Float64,1}
end
TimeData(y::Float64) =
TimeData([y])

 save_diagnostic(dname::String, d::TimeData, cname::String, c::Any, it::Integer, t::Float64) =
  pvd_add_timestep(c, field_as_vectors([0.], [0.], cname*dname, dname => (d.y), it=it, save=false), t)

mutable struct CircuitRLC
  R :: Float64
  L :: Float64
  C :: Float64
  i :: Float64
  q :: Float64
  t :: Float64
  V :: Any # in fact it is a function (Float64) ↦ (Float64)
  ext :: Any # should be some abstract type
end
CircuitRLC(i0::Number, q0::Number, t0::Number) =
	CircuitRLC(0.0, 0.0, 0.0, i0, q0, t0, nothing, nothing)

abstract type CircuitElement end
struct Resistor <: CircuitElement
	name :: String
	val  :: Float64
end
struct Inductor <: CircuitElement
	name :: String
	val  :: Float64
end
struct Capacitor <: CircuitElement
	name :: String
	val  :: Float64
end
struct VoltageSource <: CircuitElement
	name :: String
	val  :: Any
end
struct ExternalDevice <: CircuitElement
  name :: String
  val  :: Any
end

function assign!(cir::CircuitRLC, r::Resistor)  cir.R = r.val end
function assign!(cir::CircuitRLC, l::Inductor)  cir.L = l.val end
function assign!(cir::CircuitRLC, c::Capacitor) cir.C = c.val end
function assign!(cir::CircuitRLC, v::VoltageSource) cir.V = v.val end
function assign!(cir::CircuitRLC, ext::ExternalDevice) cir.ext = ext.val end
function rlc(elements)
	t0, i0, q0 = 0, 0, 0
	cir = CircuitRLC(i0, q0, t0)
	for element in elements
		println(element)
		assign!(cir, element)
	end
	cir
end
damping_factor(R::Float64,L::Float64,C::Float64) = (R/2L) * sqrt(L*C)
damping_factor(cir::CircuitRLC) = damping_factor(cir.R, cir.L, cir.C)
resonant_frequency(L::Float64,C::Float64) = 1/sqrt(L*C)/2π
resonant_frequency(cir::CircuitRLC) = resonant_frequency(cir.L, cir.C)
macro netlist(ex::Expr)
  def_netlist(ex)
end

function def_netlist(ex::Expr)
  elements = :(CircuitElement[])
  nodes = OrderedSet{String}()
  push!(nodes, "GND")
  for arg in ex.args
    parse(arg, elements, nodes)
  end
  return elements
end

function parse(ex::LineNumberNode, elements, nodes) end

function parse(ex::Expr, elements, nodes)
    if ex.head == :tuple
      name  = ex.args[1]
      node₊ = ex.args[2]
      node₋ = ex.args[3]
      rest  = ex.args[4]

      push!(nodes, string(node₊))
      push!(nodes, string(node₋))

      element = parse_circuit_element(string(name), rest)

      push!(elements.args, element)
    end
end

function parse_circuit_element(name::String, rest)
  if startswith(name, "V")
  	v = rest
  	return :(VoltageSource($(esc(name)), $(esc(v)) ))
  end

  value = rest
  if startswith(name, "R")
  	return :(Resistor($(esc(name)), $(esc(value)) ))
  elseif startswith(name, "L")
  	return :(Inductor($(esc(name)), $(esc(value)) ))
  elseif startswith(name, "C")
  	return :(Capacitor($(esc(name)), $(esc(value)) ))
  elseif startswith(name, "EXT")
    return :(ExternalDevice($(esc(name)), $(esc(value)) ))
  end
  throw("unknown element")
end

function advance_circuit!(cir::CircuitRLC, V, Δt)
  t, v = cir.t, cir.V
  i, q = cir.i, cir.q
  R, L, C = cir.R, cir.L, cir.C

  cir.i  = (L/Δt - R/2)*i + V - v(t)
  if C > 0.0
    cir.i -=  q/C
    cir.i /= (L/Δt + R/2)
    cir.q  = q + Δt*i
  else
    cir.i /= (L/Δt + R/2)
  end
  cir.t += Δt
  @diag "i" TimeData(cir.i)
  @diag "q" TimeData(cir.q)
  @diag "V" TimeData(v(t))
  @diag "Vext" TimeData(V)
end

Base.show(io :: IO, e :: Resistor) = print(io, e.name, ": ", e.val, " Ω")
Base.show(io :: IO, e :: Inductor) = print(io, e.name, ": ", e.val*1e6, " μH")
Base.show(io :: IO, e :: Capacitor)= print(io, e.name, ": ", e.val*1e9, " nF")
Base.show(io :: IO, e :: VoltageSource) = print(io, e.name, ": ", e.val, " V") 
end