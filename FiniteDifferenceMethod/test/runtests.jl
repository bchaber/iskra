using Test, Printf
import FiniteDifferenceMethod
import RegularGrid

setup_grid() = RegularGrid.create_uniform_grid(0:0.5:1, 0:0.5:2)
setup_solver(grid) = FiniteDifferenceMethod.create_poisson_solver(grid)

function parallel_capacitor()
  grid = setup_grid()
  solv = setup_solver(grid)
  bcs  = zeros(Int8, size(grid))
  bcs[ 1, 1:5] .= 1
  bcs[ 3, 1:5] .= 2
  FiniteDifferenceMethod.apply_dirichlet(solv, grid.dof[bcs .== 1], 0)
  FiniteDifferenceMethod.apply_dirichlet(solv, grid.dof[bcs .== 2], 1)
  E = FiniteDifferenceMethod.calculate_electric_field(solv, zeros(size(grid)), 0)
end

@testset "All tests" begin
  @test parallel_capacitor()[:,:,1] ≈ -1ones(3,5)
  @test parallel_capacitor()[:,:,2] ≈  zeros(3,5) atol=1e-15
end
