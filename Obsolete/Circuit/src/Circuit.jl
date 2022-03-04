module Circuit
using DataStructures
using Diagnostics

export CircuitRLC
export @netlist
export rlc
export advance_circuit!
export resonant_frequency, damping_factor

abstract type  CircuitDevice end
abstract type  CircuitElement end
mutable struct CircuitRLC
  R :: Float64
  L :: Float64
  C :: Float64
  i :: Float64
  q :: Float64
  t :: Float64
  V :: Function
  ext :: CircuitDevice
end

struct ShortedConnection <: CircuitDevice end

CircuitRLC(i0::Number, q0::Number, t0::Number) =
	CircuitRLC(0.0, 0.0, 0.0, i0, q0, t0, t -> 0., ShortedConnection())

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
	val  :: Function
end
struct ExternalDevice <: CircuitElement
  name :: String
  val  :: CircuitDevice
end

function voltage(ext::ShortedConnection) 0.0 end
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
      arg0  = ex.args[4]

      push!(nodes, string(node₊))
      push!(nodes, string(node₋))

      element = parse_circuit_element(string(name), arg0)

      push!(elements.args, element)
    end
end

function parse_circuit_element(name::String, arg0)
  value = arg0
  if startswith(name, "V")
    return :(VoltageSource($(esc(name)), $(esc(value)) ))
  elseif startswith(name, "R")
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
  vext = voltage(cir.ext)

  cir.i  = (L/Δt - R/2)*i + vext - v(t)
  if C > 0.0
    cir.i -=  q/C
    cir.i /= (L/Δt + R/2)
    cir.q  = q + Δt*i
  else L > 0.0 || R > 0.0
    cir.i /= (L/Δt + R/2)
  end
  cir.t += Δt
  @probe "Q1"   cir.q "C"
  @probe "I1"   cir.i "A"
  @probe "V1"   v(t)  "V"
  @probe "Vext" vext  "V"
end

Base.show(io :: IO, e :: Resistor) = print(io, e.name, ": ", e.val, " Ω")
Base.show(io :: IO, e :: Inductor) = print(io, e.name, ": ", e.val*1e6, " μH")
Base.show(io :: IO, e :: Capacitor)= print(io, e.name, ": ", e.val*1e9, " nF")
Base.show(io :: IO, e :: VoltageSource) = print(io, e.name, ": ", e.val, " V") 
Base.show(io :: IO, e :: ExternalDevice) = print(io, e.name, ": ", e.val, " EXT") 
end