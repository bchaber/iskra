import FiniteDifferenceMethod
import RegularGrid

using GR
h = 0.1
g = RegularGrid.create_uniform_grid(0:h:0.5, [0.0])
e = RegularGrid.create_staggered_grid(g)
nx, ny = size(g)
e["eps"] = 10ones(nx+1, ny+1, 1)
f = h^2*100ones(size(g))
f[nx] = 0; # because it has Dirichlet BC
ps = FiniteDifferenceMethod.create_generalized_poisson_solver(g,e["eps"])
FiniteDifferenceMethod.apply_dirichlet(ps, g.dof[nx,:], 0)
FiniteDifferenceMethod.apply_neumann(ps, 1)
display(ps.A), println()
display(ps.b'), println()
phi = FiniteDifferenceMethod.calculate_electric_potential(ps, -f)
plot(g.x, phi)
readline()