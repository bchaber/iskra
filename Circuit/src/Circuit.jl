module Circuit
using DataStructures

export @netlist
export rlc

mutable struct CircuitRLC
  R :: Float64
  L :: Float64
  C :: Float64
  i :: Float64
  q :: Float64
  t :: Float64
  V :: Any # in fact it is a function (Float64) ↦ (Float64)
end
CircuitRLC(i0::Number, q0::Number, t0::Number) =
	CircuitRLC(0.0, 0.0, 0.0, i0, q0, t0, nothing)

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

function assign!(cir::CircuitRLC, r::Resistor)  cir.R = r.val end
function assign!(cir::CircuitRLC, l::Inductor)  cir.L = l.val end
function assign!(cir::CircuitRLC, c::Capacitor) cir.C = c.val end
function assign!(cir::CircuitRLC, v::VoltageSource) cir.V = v.val end
function rlc(elements)
	t0, i0, q0 = 0, 0, 0
	cir = CircuitRLC(i0, q0, t0)
	for element in elements
		println(element)
		assign!(cir, element)
	end
	cir
end

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
  end
  throw("unknown element")
end

function advance!(cir::CircuitRLC, V, Δt)
	t, v = cir.t, cir.v
	i, q = cir.i, cir.q
	R, L, C = cir.R, cir.L, cir.C

	cir.i  = (L/Δt - R/2)*i + V - v(t) - q/C
	cir.i /= (L/Δt + R/2)
	cir.q  = q + Δt*i
end

Base.show(io :: IO, e :: Resistor) = print(io, e.name, ": ", e.val, " Ω")
Base.show(io :: IO, e :: Inductor) = print(io, e.name, ": ", e.val*1e9, " μH")
Base.show(io :: IO, e :: Capacitor)= print(io, e.name, ": ", e.val*1e12, " nF")
Base.show(io :: IO, e :: VoltageSource) = print(io, e.name, ": ", e.val, " V") 
end