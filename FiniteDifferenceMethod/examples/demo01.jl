import FiniteDifferenceMethod
import RegularGrid

using GR

g = RegularGrid.create_uniform_grid(0:0.1:1, 0:0.1:1)
e = RegularGrid.create_staggered_grid(g)
e["eps"] = ones(12,12,1)
e["eps"][5:8,5:8] .= 10.
g["bcs"] = zeros(11, 11, 1)
ps = FiniteDifferenceMethod.create_generalized_poisson_solver(g,e["eps"])
g["bcs"][1:11,  1] .= 1
g["bcs"][1:11, 11] .= 2
FiniteDifferenceMethod.apply_dirichlet(ps, g["bcs"] .== 1, +1)
FiniteDifferenceMethod.apply_dirichlet(ps, g["bcs"] .== 2, -1)
phi = FiniteDifferenceMethod.calculate_electric_potential(ps, zeros(size(g)))
E   = FiniteDifferenceMethod.calculate_electric_field(ps, phi)
imshow(E[:,:,1])
savefig("Ex.png")
imshow(E[:,:,2])
savefig("Ey.png")
imshow(phi)
savefig("phi.png")
