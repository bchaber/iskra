struct PoissonSolver
  A :: AbstractArray{Float64,2}
  b :: AbstractArray{Float64,1}
  εr:: AbstractArray{Float64,3}
  Δh :: Float64

  dofs :: Dict{Symbol, AbstractArray}
end
PoissonSolver(A, b, εr, Δh) =
    PoissonSolver(A, b, εr, Δh, Dict{Symbol, AbstractArray}())

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
    dofs = collect(1:nn)
    A = zeros(nn, nn)
    b = zeros(nn)
    ϕ = reshape(dofs, nx, ny)
    for j=1:ny
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
    ps = PoissonSolver(A, b, εr, Δh)
    ps.dofs[:ϕ] = ϕ
    return ps
end

function solve(ps::PoissonSolver, f)
  return ps.A\(ps.b .+ f)
end

function apply_dirichlet(ps::PoissonSolver, nodes::BitArray{3}, ϕ0)
    A, b = ps.A, ps.b
    ϕ = ps.dofs[:ϕ]
    ids = CartesianIndices(size(nodes))

    for ij in ids[nodes]
        i, j = ij.I
        A[ϕ[i,j],:]     .= 0 # clear row
        A[ϕ[i,j],ϕ[i,j]] = 1 # ϕ(i,j)
        b[ϕ[i,j]]        = ϕ0
    end
end

# Warning: it is applied at first node, 1D version only
function apply_neumann(ps::PoissonSolver, σ0)
    A, b = ps.A, ps.b
    εr = ps.εr
    Δh = ps.Δh
    ϕ = ps.dofs[:ϕ]
    # Coefficients are doubled so that the
    # the whole space charge density contributes
    # to the Boundary Condition (instead of ρ0/2)
    A[ϕ[1,1],ϕ[1,1]] = -2εr[2,2]
    A[ϕ[1,1],ϕ[2,1]] =  2εr[2,2]
    b[ϕ[1,1]]        = -2σ0*Δh
end

function calculate_electric_potential(ps::PoissonSolver, f)
    x = solve(ps, f[:])
    ϕ = x[ps.dofs[:ϕ]]
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