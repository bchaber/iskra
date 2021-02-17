import FiniteDifferenceMethod
import RegularGrids

using GR
h = 0.1
g = RegularGrids.create_uniform_grid(0:h:0.5, [0.0])
e = RegularGrids.create_staggered_grid(g)
nx, ny = size(g)
e["eps"] = 10ones(nx+1, ny+1, 1)
f = h^2*100ones(size(g))
f[nx] = 0; # because it has Dirichlet BC
g["bcs"] = zeros(nx, ny, 1)
g["bcs"][nx,:] .= 1
ps = FiniteDifferenceMethod.create_generalized_poisson_solver(g,e["eps"])
FiniteDifferenceMethod.apply_dirichlet(ps, g["bcs"] .== 1, 0)
FiniteDifferenceMethod.apply_neumann(ps, [1], 1)
display(ps.A), println()
display(ps.b'), println()
phi = FiniteDifferenceMethod.calculate_electric_potential(ps, -f)
plot(g.x, phi)
println("Press any key..."), readline()