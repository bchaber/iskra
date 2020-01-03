mutable struct PoissonSolver
  A :: AbstractArray{Float64,2}
  b :: AbstractArray{Float64,1}
  εr:: AbstractArray{Float64,3}
  Δh :: Tuple{Float64,Float64,Float64}

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
    Δx, Δy, ~ = Δh
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
    A ./= (Δx*Δx)
    ps = PoissonSolver(A, b, εr, Δh)
    ps.dofs[:ϕ] = ϕ
    return ps
end

function solve(A::Array{Float64,2}, b::Array{Float64,1})
  return A\b
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

function segment(n::Array{Int64,1})
    s, m, e = length(n) < 2 ? ([], n, []) : (n[1], n[2:end-1], n[end])
end

function apply_neumann(ps::PoissonSolver, nodes, σ0)
    A, b = ps.A, ps.b
    εr = ps.εr
    Δx, Δy, ~ = ps.Δh
    ps.dofs[:σ] = get(ps.dofs, :σ, [])
    ϕ = ps.dofs[:ϕ]
    σ = ps.dofs[:σ]
    N = maximum(ϕ)
    n = max(N, σ...) + 1
    push!(σ, n)
    A = vcat(A, zeros(1,N))
    A = hcat(A, zeros(N+1))
    b = vcat(b, [0])
    i = 1 # left edge
    nx, ny = size(ϕ)
    # Coefficients are doubled so that the
    # the whole space charge density contributes
    # to the Boundary Condition (instead of ρ0/2)
    start, middle, stop = segment(collect(nodes))
    if start == 1 push!(middle, 1); start = [] end
    if stop == nx push!(middle,nx); stop  = [] end

    for j=middle
        A[ϕ[i,j],:] .= 0
        A[ϕ[i,j],ϕ[i,  j]] -= 2εr[i,  j]/Δx^2
        A[ϕ[i,j],ϕ[i+1,j]] += 2εr[i+1,j]/Δx^2
        A[ϕ[i,j],σ[end]]   += 2/Δx
    end
    
    for j=start
        A[ϕ[i,j],:] .= 0
        A[ϕ[i,j],ϕ[i,j]]   -= 4εr[i,j]/3Δx^2
        A[ϕ[i,j],ϕ[i,j]]   -= 4εr[i,j]/3Δy^2
        A[ϕ[i,j],ϕ[i+1,j]] += 4εr[i,j]/3Δx^2
        A[ϕ[i,j],ϕ[i,j-1]] += 4εr[i,j]/3Δy^2
        A[ϕ[i,j],σ[end]]   += (2/3)*(Δx+Δy)/(Δx*Δy)
    end

    for j=stop
        A[ϕ[i,j],:] .= 0
        A[ϕ[i,j],ϕ[i,j]]   -= 4εr[i,j]/3Δx^2
        A[ϕ[i,j],ϕ[i,j]]   -= 4εr[i,j]/3Δy^2
        A[ϕ[i,j],ϕ[i+1,j]] += 4εr[i,j]/3Δx^2
        A[ϕ[i,j],ϕ[i,j+1]] += 4εr[i,j]/3Δy^2
        A[ϕ[i,j],σ[end]]   += (2/3)*(Δx+Δy)/(Δx*Δy)
    end

    A[σ[end],σ[end]] = 1
    b[σ[end]]        = σ0
    ps.A, ps.b = A, b
end

function get_rhs(ps::PoissonSolver, s::Symbol, dof::Int64) =
    view(ps.b, ps.dofs[s][dof])

function calculate_electric_potential(ps::PoissonSolver, f)
    A, b = ps.A, ps.b
    dofs = ps.dofs[:ϕ]
    b[dofs[2:end-1,2:end-1]] .= f[dofs[2:end-1,2:end-1]]
    x = solve(A, b)
    ϕ = x[ps.dofs[:ϕ]]
    return ϕ
end

function calculate_electric_field(ps::PoissonSolver, ϕ)
    nx, ny = size(ϕ)
    calculate_electric_field!(ps, ϕ, zeros(nx, ny, 3))
end

function calculate_electric_field!(ps::PoissonSolver, ϕ, E)
    nx, ny = size(ϕ)
    Δx, Δy, ~ = ps.Δh

    E[2:nx-1,:,1] = (ϕ[1:nx-2,:] - ϕ[3:nx,:])/2Δx # central difference on internal nodes
    E[:,2:ny-1,2] = (ϕ[:,1:ny-2] - ϕ[:,3:ny])/2Δy # central difference on internal nodes
    E[1 ,:,1]  = (ϕ[1,:]    - ϕ[ 2,:])/Δx  #  forward difference on x=0
    E[nx,:,1]  = (ϕ[nx-1,:] - ϕ[nx,:])/Δx  # backward difference on x=Lx
    E[:, 1,2]  = (ϕ[:,1]    - ϕ[:, 2])/Δy  #  forward difference on y=0
    E[:,ny,2]  = (ϕ[:,ny-1] - ϕ[:,ny])/Δy  # backward difference on y=Ly

    return E
end