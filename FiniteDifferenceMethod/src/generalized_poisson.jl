abstract type FieldBoundary end
struct Open <: FieldBoundary end
struct Other <: FieldBoundary end
struct Periodic <: FieldBoundary end
struct Reflecting <: FieldBoundary end
struct Top end
struct Right end
struct Bottom end
struct Left end

mutable struct PoissonSolver{CS, D}
  A :: AbstractArray{Float64,2}
  b :: AbstractArray{Float64,1}
  x :: AbstractArray{Float64,1}
  εr:: AbstractArray{Float64,D}
  ε0:: Float64
  Δh :: NTuple{D, Float64}
  bnds :: Dict{Symbol, FieldBoundary}
  dofs :: Dict{Symbol, AbstractArray}
end

function create_poisson_solver(grid::CartesianGrid{2}, ε0::Float64)
    nx, ny = size(grid)
    create_generalized_poisson_solver(grid, ones(nx+1, ny+1), ε0)
end

# Implement generalized Poisson solver from (eq. 37-40):
# https://my.ece.utah.edu/~ece6340/LECTURES/Feb1/Nagel 2012 - Solving the Generalized Poisson Equation using FDM.pdf
function create_generalized_poisson_solver(grid::CartesianGrid{2}, εr::Array{Float64,2}, ε0::Float64)
    nx, ny = size(grid)
    nn = nx⋅ny
    Δx, Δy = grid.Δh
    A = zeros(nn, nn)
    b = zeros(nn)
    x = zeros(nn)
    ϕ = reshape(collect(1:nn), nx, ny)
    ρ = reshape(collect(1:nn), nn)
    for j=1:ny
        for i=1:nx
            if i < nx
                A[ϕ[i,j],ϕ[i,j]]   -= 0.5εr[i+1,j+1] + 0.5εr[i+1,j]   # ϕ(i,j)
                A[ϕ[i,j],ϕ[i+1,j]] += 0.5εr[i+1,j+1] + 0.5εr[i+1,j]   # ϕ(i+1,j)
            end
            if i > 1
                A[ϕ[i,j],ϕ[i,j]]   -= 0.5εr[i,j]     + 0.5εr[i,j+1]   # ϕ(i,j)
                A[ϕ[i,j],ϕ[i-1,j]] += 0.5εr[i,j]     + 0.5εr[i,j+1]   # ϕ(i-1,j)
            end
            if j < ny
                A[ϕ[i,j],ϕ[i,j]]   -= 0.5εr[i+1,j+1] + 0.5εr[i,j+1]   # ϕ(i,j)
                A[ϕ[i,j],ϕ[i,j+1]] += 0.5εr[i+1,j+1] + 0.5εr[i,j+1]   # ϕ(i,j+1)
            end
            if j > 1
                A[ϕ[i,j],ϕ[i,j]]   -= 0.5εr[i+1,j]   + 0.5εr[i,j]     # ϕ(i,j)
                A[ϕ[i,j],ϕ[i,j-1]] += 0.5εr[i+1,j]   + 0.5εr[i,j]     # ϕ(i,j-1)
            end
        end
    end
    A ./= (Δx*Δx) # TODO: it will break when Δx ≠ Δy
    dofs = Dict{Symbol, AbstractArray}(:ϕ => ϕ, :ρ => ρ)
    bnds = Dict{Symbol, FieldBoundary}(:left => Open(), :right => Open(),
                                       :bottom => Open(), :top => Open())
    ps = PoissonSolver{:xy, 2}(A, b, x, εr, ε0, (Δx, Δy), bnds, dofs)
    return ps
end

function solve(A, b)
  return A\b
end

function apply_dirichlet(ps::PoissonSolver{:xy, 2}, nodes::BitArray{2}, ϕ0)
    A, b = ps.A, ps.b
    ϕ, ρ = ps.dofs[:ϕ], ps.dofs[:ρ]
    ids = CartesianIndices(size(nodes))
    if ϕ0 ≠ 0
        println("ERROR: non-zero Dirichlet Boundary Condition")
    end
    for ij in ids[nodes]
        i, j = ij.I
        A[ϕ[i,j],:]     .= 0 # clear row
        A[ϕ[i,j],ϕ[i,j]] = 1 # ϕ(i,j)
        b[ϕ[i,j]]        = ϕ0
        setdiff!(ρ, ϕ[i,j])
    end
end

function add_new_dof(ps::PoissonSolver, symbol::Symbol)
    N = length(ps.b)
    ps.dofs[symbol] = s = get(ps.dofs, symbol, [])
    push!(s,  N+1)
    
    ps.A = vcat(ps.A, zeros(1,N))
    ps.A = hcat(ps.A, zeros(N+1))
    ps.b = vcat(ps.b, [0.])
    
    ps.x = vcat(ps.x, [0.])
    ps.A[s[end],s[end]] = 1
    ps.b[s[end]]        = 0
    return length(s)
end

# Coefficients are doubled so that the
# the whole space charge density contributes
# to the Boundary Condition (instead of ρ0/2)
function apply_neumann(ps::PoissonSolver{:xy, 2}, nodes, dof)
    εr = ps.εr
    Δx, Δy = ps.Δh
    ϕ = ps.dofs[:ϕ]
    σ = ps.dofs[:σ]
    nx, ny = size(ϕ)
    A, b = ps.A, ps.b

    for c in findall(nodes)
        i, j = c.I
        a = (i >  1) ? nodes[i-1,j] : false
        b = (i < nx) ? nodes[i+1,j] : false
        c = (j >  1) ? nodes[i,j-1] : true
        d = (j < ny) ? nodes[i,j+1] : true

        if c == d == true &&
           a == b == false
            i′ = (i == 1) ?   i+1 : i-1
            A[ϕ[i,j],:] .= 0
            A[ϕ[i,j],ϕ[i, j]] -= 2εr[i, j]/Δx
            A[ϕ[i,j],ϕ[i′,j]] += 2εr[i′,j]/Δx
            A[ϕ[i,j],σ[dof]]  += 2
        end
        if c != d
            i′ = (i == 1) ? i+1 : i-1
            j′ = (c) ?      j+1 : j-1
            A[ϕ[i,j],:] .= 0
            A[ϕ[i,j],ϕ[i, j]]  -= 4εr[i,j]/3Δx^2
            A[ϕ[i,j],ϕ[i, j]]  -= 4εr[i,j]/3Δy^2
            A[ϕ[i,j],ϕ[i′,j]]  += 4εr[i,j]/3Δx^2
            A[ϕ[i,j],ϕ[i,j′]]  += 4εr[i,j]/3Δy^2
            A[ϕ[i,j],σ[dof]]   += (2/3)*(Δx+Δy)/(Δx*Δy)
        end
    end
end

function apply_periodic(ps::PoissonSolver{:xy, 1})
    εr,A = ps.εr, ps.A
    Δx,  = ps.Δh
    ϕ, ρ = ps.dofs[:ϕ], ps.dofs[:ρ], ps.A
    nx,  = size(ϕ)
    ps.bnds[:left] = Periodic()
    ps.bnds[:right] = Periodic()
    
    A[ϕ[1], ϕ[ 1]] -= (0.5εr[2] + 0.5εr[1])/Δx^2
    A[ϕ[1], ϕ[nx]] += (0.5εr[2] + 0.5εr[1])/Δx^2

    A[ϕ[nx],ϕ[nx]] -= (0.5εr[2] + 0.5εr[1])/Δx^2
    A[ϕ[nx],ϕ[ 1]] += (0.5εr[2] + 0.5εr[1])/Δx^2
end

function apply_periodic(ps::PoissonSolver{:xy, 2}, axis)
    εr,A = ps.εr, ps.A
    Δx, Δy = ps.Δh
    ϕ, ρ = ps.dofs[:ϕ], ps.dofs[:ρ], ps.A
    nx, ny = size(ϕ)
    if axis == 1
        ps.bnds[:left] = Periodic()
        ps.bnds[:right] = Periodic()
        for j in (1, ny)
            for i=1:nx
                if j == ny
                    A[ϕ[i,j],ϕ[i,j]] -= (0.5εr[i+1,j+1] + 0.5εr[i,j+1])/Δx^2
                    A[ϕ[i,j],ϕ[i,1]] += (0.5εr[i+1,j+1] + 0.5εr[i,j+1])/Δx^2
                end
                if j == 1
                    A[ϕ[i,j],ϕ[i, j]] -= (0.5εr[i+1,j]   + 0.5εr[i,j])/Δx^2
                    A[ϕ[i,j],ϕ[i,ny]] += (0.5εr[i+1,j]   + 0.5εr[i,j])/Δx^2
                end
            end
        end
    end

    if axis == 2
        ps.bnds[:top] = Periodic()
        ps.bnds[:bottom] = Periodic()
        for j=1:ny
            for i in (1, nx)
                if i == nx
                    A[ϕ[i,j],ϕ[i,j]] -= (0.5εr[i+1,j+1] + 0.5εr[i+1,j])/Δy^2
                    A[ϕ[i,j],ϕ[1,j]] += (0.5εr[i+1,j+1] + 0.5εr[i+1,j])/Δy^2
                end
                if i == 1
                    A[ϕ[i,j],ϕ[i, j]] -= (0.5εr[i,j]     + 0.5εr[i,j+1])/Δy^2
                    A[ϕ[i,j],ϕ[nx,j]] += (0.5εr[i,j]     + 0.5εr[i,j+1])/Δy^2
                end
            end
        end
    end
end

get_solution(ps::PoissonSolver, s::Symbol, i::Int64, j::Int64) =
    view(ps.x, ps.dofs[s][i,j])
get_solution(ps::PoissonSolver, s::Symbol, i::Int64) =
    view(ps.x, ps.dofs[s][i])
get_rhs(ps::PoissonSolver, s::Symbol, i::Int64, j::Int64) =
    view(ps.b, ps.dofs[s][i,j])
get_rhs(ps::PoissonSolver, s::Symbol, i::Int64) =
    view(ps.b, ps.dofs[s][i])


update!(::Open,     ::Top,    f, ϕ, nx, ny) = nothing #f[ϕ[nx,:]] .= 2.0f[ϕ[nx,:]]
update!(::Open,     ::Right,  f, ϕ, nx, ny) = nothing #f[ϕ[:,ny]] .= 2.0f[ϕ[:,ny]]
update!(::Open,     ::Bottom, f, ϕ, nx, ny) = nothing #f[ϕ[1, :]] .= 2.0f[ϕ[1, :]]
update!(::Open,     ::Left,   f, ϕ, nx, ny) = nothing #f[ϕ[:, 1]] .= 2.0f[ϕ[:, 1]]

update!(::Periodic, ::Top,    f, ϕ, nx, ny) = f[ϕ[nx,:]] .= f[ϕ[1, :]] .+ f[ϕ[nx, :]]
update!(::Periodic, ::Right,  f, ϕ, nx, ny) = f[ϕ[:,ny]] .= f[ϕ[:, 1]] .+ f[ϕ[:, ny]]
update!(::Periodic, ::Bottom, f, ϕ, nx, ny) = f[ϕ[1, :]] .= f[ϕ[1, :]] .+ f[ϕ[nx, :]]
update!(::Periodic, ::Left,   f, ϕ, nx, ny) = f[ϕ[:, 1]] .= f[ϕ[:, 1]] .+ f[ϕ[:, ny]]

function update_source_term!(ps::PoissonSolver{:xy, 2}, f, ϕ)
    nx, ny = size(ϕ)
    update!(ps.bnds[:top],    Top(),    f, ϕ, nx, ny)
    update!(ps.bnds[:right],  Right(),  f, ϕ, nx, ny)
    update!(ps.bnds[:bottom], Bottom(), f, ϕ, nx, ny)
    update!(ps.bnds[:left],   Left(),   f, ϕ, nx, ny)
end

function calculate_electric_potential(ps::PoissonSolver{CS, 2}, f) where CS
    A, b, x, ε0 = ps.A, ps.b, ps.x, ps.ε0
    ϕ, ρ = ps.dofs[:ϕ], ps.dofs[:ρ]
    update_source_term!(ps, f, ϕ)
    b[ρ] .= f[ρ]
    x .= solve(A, b/ε0)
    return x[ϕ] # we could use (or even reuse!) a view here
end

function calculate_electric_field(ps::PoissonSolver{CS,2}, ϕ) where CS
    nx, ny = size(ϕ)
    calculate_electric_field!(ps, ϕ, zeros(nx, ny, 3))
end

function calculate_electric_field!(ps::PoissonSolver{:xy, 2}, ϕ, E)
    nx, ny = size(ϕ)
    Δx, Δy = ps.Δh

    E[2:nx-1,:,1] = (ϕ[1:nx-2,:] - ϕ[3:nx,:])/2Δx # central difference on internal nodes
    E[:,2:ny-1,2] = (ϕ[:,1:ny-2] - ϕ[:,3:ny])/2Δy # central difference on internal nodes
    E[1 ,:,1]  = (ϕ[1,:]    - ϕ[ 2,:])/Δx  #  forward difference on x=0
    E[nx,:,1]  = (ϕ[nx-1,:] - ϕ[nx,:])/Δx  # backward difference on x=Lx
    E[:, 1,2]  = (ϕ[:,1]    - ϕ[:, 2])/Δy  #  forward difference on y=0
    E[:,ny,2]  = (ϕ[:,ny-1] - ϕ[:,ny])/Δy  # backward difference on y=Ly

    return E
end