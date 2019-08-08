module FiniteDifferenceMethod
    export calculate_electric_field

    using RegularGrid
    using LinearAlgebra
    using Printf:@sprintf,println

    mutable struct LinearSystemOfEquations
      A :: AbstractArray{Float64,2}
      b :: AbstractArray{Float64,1}
    end
 
    struct PoissonSolver
      Δh :: Float64
    end
   
    lse = LinearSystemOfEquations(zeros(0,0), zeros(0))

import Diagnostics
struct GridData <: Diagnostics.DiagnosticData
  u :: Array{Float64,2}
 sp :: Array{Float64,1}
 or :: Array{Float64,1}
 it :: Integer
end

import PlotVTK: pvd_add_timestep, heatmap
Diagnostics.save_diagnostic(dname::String, d::GridData, cname::String, c::Any, it::Integer) =
  pvd_add_timestep(c, heatmap(dname => d.u, dname, spacing=d.sp, origin=d.or, it=it, save=false), it)

function create_poisson_solver(grid::UniformGrid)
    nx, ny = size(grid)
    nn = nx⋅ny
    Δh = grid.Δh

    A = lse.A = zeros(nn, nn)
    b = lse.b = zeros(nn)
    ϕ = reshape(1:nn, nx, ny)
    # set regular stencil on internal nodes
    for j=2:ny-1             # only internal nodes
        for i=2:nx-1
            A[ϕ[i,j],ϕ[i,j]]   +=-2/Δh^2 # ϕ(i,j)
            A[ϕ[i,j],ϕ[i-1,j]] += 1/Δh^2 # ϕ(i-1,j)
            A[ϕ[i,j],ϕ[i+1,j]] += 1/Δh^2 # ϕ(i+1,j)

            A[ϕ[i,j],ϕ[i,j]]   +=-2/Δh^2 # ϕ(i,j)
            A[ϕ[i,j],ϕ[i,j-1]] += 1/Δh^2 # ϕ(i,j-1)
            A[ϕ[i,j],ϕ[i,j+1]] += 1/Δh^2 # ϕ(i,j+1)
        end
    end

    # y=0
    for j=1
        for i=2:nx-1
            A[ϕ[i,j],ϕ[i,j]]   +=-2/Δh^2 # ϕ(i,j)
            A[ϕ[i,j],ϕ[i-1,j]] += 1/Δh^2 # ϕ(i-1,j)
            A[ϕ[i,j],ϕ[i+1,j]] += 1/Δh^2 # ϕ(i+1,j)

            A[ϕ[i,j],ϕ[i,j]]   +=-1/Δh^2 # ϕ(i,j)
            A[ϕ[i,j],ϕ[i,j+1]] += 1/Δh^2 # ϕ(i,j+1)
        end
    end

    # y=Ly
    for j=ny
        for i=2:nx-1
            A[ϕ[i,j],ϕ[i,j]]   +=-2/Δh^2 # ϕ(i,j)
            A[ϕ[i,j],ϕ[i-1,j]] += 1/Δh^2 # ϕ(i-1,j)
            A[ϕ[i,j],ϕ[i+1,j]] += 1/Δh^2 # ϕ(i+1,j)

            A[ϕ[i,j],ϕ[i,j]]   +=-1/Δh^2 # ϕ(i,j)
            A[ϕ[i,j],ϕ[i,j-1]] += 1/Δh^2 # ϕ(i,j-1)
        end
    end

    # x=Lx
    for i=nx
        for j=2:ny-1
            A[ϕ[i,j],ϕ[i,j]]   +=-1/Δh^2 # ϕ(i,j)
            A[ϕ[i,j],ϕ[i-1,j]] += 1/Δh^2 # ϕ(i-1,j)

            A[ϕ[i,j],ϕ[i,j]]   +=-2/Δh^2 # ϕ(i,j)
            A[ϕ[i,j],ϕ[i,j-1]] += 1/Δh^2 # ϕ(i,j-1)
            A[ϕ[i,j],ϕ[i,j+1]] += 1/Δh^2 # ϕ(i,j+1)
        end
    end

    # x=0
    for i=1
        for j=2:ny-1
            A[ϕ[i,j],ϕ[i,j]]   +=-1/Δh^2 # ϕ(i,j)
            A[ϕ[i,j],ϕ[i+1,j]] += 1/Δh^2 # ϕ(i+1,j)

            A[ϕ[i,j],ϕ[i,j]]   +=-2/Δh^2 # ϕ(i,j)
            A[ϕ[i,j],ϕ[i,j-1]] += 1/Δh^2 # ϕ(i,j-1)
            A[ϕ[i,j],ϕ[i,j+1]] += 1/Δh^2 # ϕ(i,j+1)
        end
    end
    
    for i=1
        for j=1
            A[ϕ[i,j],ϕ[i,j]]   +=-1/Δh^2 # ϕ(i,j)
            A[ϕ[i,j],ϕ[i+1,j]] += 1/Δh^2 # ϕ(i+1,j)

            A[ϕ[i,j],ϕ[i,j]]   +=-1/Δh^2 # ϕ(i,j)
            A[ϕ[i,j],ϕ[i,j+1]] += 1/Δh^2 # ϕ(i,j+1)
        end
    end
    
    for i=1
        for j=ny
            A[ϕ[i,j],ϕ[i,j]]   +=-1/Δh^2 # ϕ(i,j)
            A[ϕ[i,j],ϕ[i+1,j]] += 1/Δh^2 # ϕ(i+1,j)

            A[ϕ[i,j],ϕ[i,j]]   +=-1/Δh^2 # ϕ(i,j)
            A[ϕ[i,j],ϕ[i,j-1]] += 1/Δh^2 # ϕ(i,j-1)
        end
    end
    
    for i=nx
        for j=1
            A[ϕ[i,j],ϕ[i,j]]   +=-1/Δh^2 # ϕ(i,j)
            A[ϕ[i,j],ϕ[i-1,j]] += 1/Δh^2 # ϕ(i-1,j)

            A[ϕ[i,j],ϕ[i,j]]   +=-1/Δh^2 # ϕ(i,j)
            A[ϕ[i,j],ϕ[i,j+1]] += 1/Δh^2 # ϕ(i,j+1)
        end
    end
    
    for i=nx
        for j=ny
            A[ϕ[i,j],ϕ[i,j]]   +=-1/Δh^2 # ϕ(i,j)
            A[ϕ[i,j],ϕ[i-1,j]] += 1/Δh^2 # ϕ(i-1,j)

            A[ϕ[i,j],ϕ[i,j]]   +=-1/Δh^2 # ϕ(i,j)
            A[ϕ[i,j],ϕ[i,j-1]] += 1/Δh^2 # ϕ(i,j-1)
        end
    end
    return PoissonSolver(Δh)
end

function solve!()
  return lse.A\lse.b
end

function apply_dirichlet(::PoissonSolver, dofs, value)
    A, b = lse.A, lse.b
    for dof=dofs
        ϕ = reshape(1:length(A), size(A))
        i, j = findfirst(x -> x == dof, ϕ).I
        A[ϕ[i,j],:]     .= 0 # clear row
        A[ϕ[i,j],ϕ[i,j]] = 1 # ϕ(i,j)
        b[ϕ[i,j]]        = value
    end
end

function calculate_electric_potential(ρ)
    x = solve!()
    ϕ = reshape(x, size(ρ))
    return ϕ
end

function calculate_electric_field(ps::PoissonSolver, ρ, it)
    ϕ = calculate_electric_potential(ρ)
    nx, ny = size(ϕ)
    Ex = zeros(nx, ny)
    Ey = zeros(nx, ny)
    Δh = ps.Δh

    Ex[2:nx-1,:] = ϕ[1:nx-2,:] - ϕ[3:nx,:]  # central difference on internal nodes
    Ey[:,2:ny-1] = ϕ[:,1:ny-2] - ϕ[:,3:ny]  # central difference on internal nodes
    Ex[1 ,:]  = 2*(ϕ[1,:]    - ϕ[ 2,:])     #  forward difference on x=0
    Ex[nx,:]  = 2*(ϕ[nx-1,:] - ϕ[nx,:])     # backward difference on x=Lx
    Ey[:, 1]  = 2*(ϕ[:,1]    - ϕ[:, 2])     #  forward difference on y=0
    Ey[:,ny]  = 2*(ϕ[:,ny-1] - ϕ[:,ny])     # backward difference on y=Ly

    Ex = Ex/2Δh
    Ey = Ey/2Δh

    origin, spacing = [Δh,Δh], [0,0]
    Diagnostics.register_diagnostic("ϕ",  GridData(ϕ,  origin, spacing, it))
    Diagnostics.register_diagnostic("Ex", GridData(Ex, origin, spacing, it))
    Diagnostics.register_diagnostic("Ey", GridData(Ey, origin, spacing, it))

    return Ex, Ey
end
end
