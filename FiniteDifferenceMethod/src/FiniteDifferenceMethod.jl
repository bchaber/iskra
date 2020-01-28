module FiniteDifferenceMethod
	export create_poisson_solver
    export calculate_electric_field
    export calculate_electric_field!
    export calculate_electric_potential
    export calculate_advection_diffusion
    using RegularGrid
    using LinearAlgebra
 
    include("generalized_poisson.jl")
    include("advection_diffusion.jl")
end
