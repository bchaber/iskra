module FiniteDifferenceMethod
    export calculate_electric_field
    export calculate_electric_field!
    export calculate_electric_potential
    export calculate_advection_diffusion
    using RegularGrid
    using LinearAlgebra

    mutable struct LinearSystemOfEquations
      A :: AbstractArray{Float64,2}
      b :: AbstractArray{Float64,1}
    end
 
    struct PoissonSolver
      Δh :: Float64
    end
   
    lse = LinearSystemOfEquations(zeros(0,0), zeros(0))

function create_poisson_solver(grid::UniformGrid)
    nx, ny = size(grid)
    create_generalized_poisson_solver(grid, ones(nx+1, ny+1, 1))
end

# Implement generalized Poisson solver from (eq. 37-40):
# https://my.ece.utah.edu/~ece6340/LECTURES/Feb1/Nagel 2012 - Solving the Generalized Poisson Equation using FDM.pdf
function create_generalized_poisson_solver(grid::UniformGrid, εr::Array{Float64,3})
    nx, ny = size(grid)
    nn = nx⋅ny
    Δh = grid.Δh

    A = lse.A = zeros(nn, nn)
    b = lse.b = zeros(nn)
    ϕ = reshape(1:nn, nx, ny)
    # set regular stencil on internal nodes
    for j=1:ny             # only internal nodes
        for i=1:nx
            if i < nx
                A[ϕ[i,j],ϕ[i,j]]   += -0.5εr[i+1,j+1] - 0.5εr[i+1,j]   # ϕ(i,j)
                A[ϕ[i,j],ϕ[i+1,j]] +=  0.5εr[i+1,j+1] + 0.5εr[i+1,j]   # ϕ(i+1,j)
            end
            if i > 1
                A[ϕ[i,j],ϕ[i,j]]   += -0.5εr[i,j]   - 0.5εr[i,j+1]     # ϕ(i,j)
                A[ϕ[i,j],ϕ[i-1,j]] +=  0.5εr[i,j]   + 0.5εr[i,j+1]     # ϕ(i-1,j)
            end
            if j < ny
                A[ϕ[i,j],ϕ[i,j]]   += -0.5εr[i+1,j+1] - 0.5εr[i,j+1]   # ϕ(i,j)
                A[ϕ[i,j],ϕ[i,j+1]] +=  0.5εr[i+1,j+1] + 0.5εr[i,j+1]   # ϕ(i,j+1)
            end
            if j > 1
                A[ϕ[i,j],ϕ[i,j]]   += -0.5εr[i+1,j]   - 0.5εr[i,j]     # ϕ(i,j)
                A[ϕ[i,j],ϕ[i,j-1]] +=  0.5εr[i+1,j]   + 0.5εr[i,j]     # ϕ(i,j-1)
            end
        end
    end
    return PoissonSolver(Δh)
end

function solve(f)
  return lse.A\(lse.b .+ f)
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

function calculate_advection_diffusion(n, v, Δh, Δt)
    nx, ny = size(n)
    vx = view(v,:,:,1)
    vy = view(v,:,:,2)
    Δn = zeros(nx, ny)

    for i=2:nx-1
      for j=2:ny-1
        # diffusion
        Δn[i,j] += Δt*1*(n[i-1,j] + n[i,j-1] - 4n[i,j] + n[i+1,j] + n[i,j+1])/Δh^2
        # advection
        Δn[i,j] -= Δt*n[i,j] * (vx[i+1,j] - vx[i-1,j] + vy[i,j+1] - vy[i,j-1])/2Δh
        Δn[i,j] -= Δt*(vx[i,j] * (n[i+1,j] - n[i-1,j]) + vy[i,j] * (n[i,j+1] - n[i,j-1]))/Δh 
      end
    end

    return Δn
end

function calculate_electric_potential(ps::PoissonSolver, f)
    x = solve(f[:])
    ϕ = reshape(x, size(f))
    return ϕ
end

function calculate_electric_field(ps::PoissonSolver, ϕ)
    nx, ny = size(ϕ)
    calculate_electric_field!(ps, ϕ, zeros(nx, ny, 3))
end

function calculate_electric_field!(ps::PoissonSolver, ϕ, E)
    nx, ny = size(ϕ)
    Δh = ps.Δh

    E[2:nx-1,:,1] = ϕ[1:nx-2,:] - ϕ[3:nx,:]  # central difference on internal nodes
    E[:,2:ny-1,2] = ϕ[:,1:ny-2] - ϕ[:,3:ny]  # central difference on internal nodes
    E[1 ,:,1]  = 2*(ϕ[1,:]    - ϕ[ 2,:])     #  forward difference on x=0
    E[nx,:,1]  = 2*(ϕ[nx-1,:] - ϕ[nx,:])     # backward difference on x=Lx
    E[:, 1,2]  = 2*(ϕ[:,1]    - ϕ[:, 2])     #  forward difference on y=0
    E[:,ny,2]  = 2*(ϕ[:,ny-1] - ϕ[:,ny])     # backward difference on y=Ly

    E ./=  2Δh

    return E
end
end
