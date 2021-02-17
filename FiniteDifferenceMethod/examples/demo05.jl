import FiniteDifferenceMethod
import Circuit: @netlist, advance_circuit!, rlc
import RegularGrids

σ0, ρ0 = 0.0, 0.0
m = 81
h = 0.5/(m-1)
g = RegularGrids.create_uniform_grid(0:h:(m-1)*h, [0.0])
e = RegularGrids.create_staggered_grid(g)
nx, ny = size(g)
se = round(Int64, nx/3)+1
ee = round(Int64,2nx/3)-1
e["eps"] = 10ones(nx+1, ny+1, 1)
ρ = zeros(size(g))
ρ[se:ee] .= ρ0
ρ[nx] = 0; # because it has Dirichlet BC
g["bcs"] = zeros(nx, ny, 1)
g["bcs"][nx,:] .= 1
# Field solver
ps = FiniteDifferenceMethod.create_generalized_poisson_solver(g,e["eps"])
FiniteDifferenceMethod.apply_dirichlet(ps, g["bcs"] .== 1, 0)
FiniteDifferenceMethod.apply_neumann(ps, [1], σ0)
ϕ = FiniteDifferenceMethod.calculate_electric_potential(ps, -ρ)
# Circuit solver
v(t) = 5e-9 < t < 10e-9 ? 1. : 0.
cir = rlc(@netlist begin
	V1, VCC, GND, v
	L1, NOD, VCC, 5e-9
	C1, NOD, VCC, 1e-9
	R1, GND, NOD, 4.5
end)
ζ = (cir.R/2cir.L) * sqrt(cir.L*cir.C)
println("damping factor ζ = ", ζ)
# Field-Circuit integration
data = zeros(1000, 4)
Δt = 0.1e-9
Δx, Δy, Δz = g.Δh
using GR

for i=1:size(data,1)
	global ϕ
	data[i,:] = [cir.t, cir.i, cir.q, ϕ[1]]
	V = ϕ[1] - ϕ[nx]
	advance_circuit!(cir, V, Δt)
	dσ = -Δt*cir.i/(Δy*Δz)
	ps.b[ps.dofs[:σ][1]] += dσ
	ϕ = FiniteDifferenceMethod.calculate_electric_potential(ps, -ρ)
end
t, ϕ0 = data[:,1], data[:,4]
plot(t, ϕ0, t, maximum(ϕ0)*v.(t))
print("press any key..."), readline()