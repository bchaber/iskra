import FiniteDifferenceMethod
import Circuit: @netlist, advance!, rlc
import RegularGrid

σ = 0.0
h = 0.5/15
g = RegularGrid.create_uniform_grid(0:h:15h, [0.0])
e = RegularGrid.create_staggered_grid(g)
nx, ny = size(g)
e["eps"] = 10ones(nx+1, ny+1, 1)
ρ = 100ones(size(g))
ρ[nx] = 0; # because it has Dirichlet BC
g["bcs"] = zeros(nx, ny, 1)
g["bcs"][nx,:] .= 1
# Field solver
ps = FiniteDifferenceMethod.create_generalized_poisson_solver(g,e["eps"])
FiniteDifferenceMethod.apply_dirichlet(ps, g["bcs"] .== 1, 0)
FiniteDifferenceMethod.apply_neumann(ps, σ)
ϕ = FiniteDifferenceMethod.calculate_electric_potential(ps, -ρ*h^2)
# Circuit solver
v(t) = sin(2π*200*t)
cir = rlc(@netlist begin
	V1, VCC, GND, t -> 1
	L1, NOD, VCC, 1e-3
	C1, NOD, VCC, 50e-12
	R1, GND, NOD, 200
end)
# Field-Circuit integration
data = zeros(10000, 3)
Δt = 0.01e-9
Δx, Δy, Δz = g.Δh
using GR

for i=1:size(data,1)
	global ϕ
	data[i,:] = [cir.t, cir.i, cir.q]
	V = 0.0
	advance!(cir, V, Δt)
	dσ = -Δt*cir.i/(Δy*Δz)
	ps.b[ps.dofs[:σ][1]] += dσ
	ϕ = FiniteDifferenceMethod.calculate_electric_potential(ps, -ρ*h^2)
end
plot(data[:,1], data[:,3])
print("press any key..."), readline()