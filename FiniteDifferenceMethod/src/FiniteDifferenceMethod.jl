module FiniteDifferenceMethod
    export setup, solve

    using RegularGrid
    using LinearAlgebra
    using Printf:@sprintf,println

    function setup(grid::UniformGrid)
        global A
        global b

        nx, ny = size(grid)
        nn = nx⋅ny
        Δh = grid.Δh

        A = zeros(nn, nn)
        b = zeros(nn, 1)
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
    end

    function solve!()
      return A\b
    end

    function apply_dirichlet(dofs, value)
        global A
        global b

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
    
    function calculate_electric_field(ϕ, Δh)
        nx, ny = size(ϕ)
        Ex = zeros(nx, ny)
        Ey = zeros(nx, ny)

        Ex[2:nx-1,:] = ϕ[1:nx-2,:] - ϕ[3:nx,:]  # central difference on internal nodes
        Ey[:,2:ny-1] = ϕ[:,1:ny-2] - ϕ[:,3:ny]  # central difference on internal nodes
        Ex[1 ,:]  = 2*(ϕ[1,:]    - ϕ[ 2,:])     #  forward difference on x=0
        Ex[nx,:]  = 2*(ϕ[nx-1,:] - ϕ[nx,:])     # backward difference on x=Lx
        Ey[:, 1]  = 2*(ϕ[:,1]    - ϕ[:, 2])     #  forward difference on y=0
        Ey[:,ny]  = 2*(ϕ[:,ny-1] - ϕ[:,ny])     # backward difference on y=Ly
    
        Ex = Ex/2Δh
        Ey = Ey/2Δh
    
        return Ex, Ey
    end
end
