module FiniteDifferenceMethod
    export apply_dirichlet, apply_neumann, apply_periodic
    export get_rhs, get_solution
    export create_poisson_solver
    export create_generalized_poisson_solver
    export calculate_electric_field
    export calculate_electric_field!
    export calculate_electric_potential
    export calculate_advection_diffusion
    export calculate_magnetic_field
    using RegularGrids
    using LinearAlgebra
 
    include("generalized_poisson.jl")
    include("advection_diffusion.jl")
end
