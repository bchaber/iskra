module RegularGrid
    export create_uniform_grid, create_staggered_grid
    export create_axial_grid
    export UniformGrid, CartesianGrid, AxialGrid
    export cell_volume

    struct UniformGrid{C, D}
      data::Dict{String,AbstractArray}
         n::NTuple{D,Integer}
        Δh::NTuple{D,Float64}
      node::Array{Float64,D}
    coords::NTuple{D,Array{Float64,D}}
    origin::NTuple{D,Float64}
    end

    CartesianGrid{D} = UniformGrid{:xy, D}
    AxialGrid{D}     = UniformGrid{:rz, D}

    Base.size(g::UniformGrid) = g.n
    Base.length(g::UniformGrid) = prod(g.n)
    Base.getindex(g::UniformGrid, k) = g.data[k]
    Base.setindex!(g::UniformGrid, v, k) =
        length(v) ≠ length(g) ? error("Size mismatch!") : setindex!(g.data, v, k)

    function cell_volume(g::CartesianGrid{2})
        Δx, Δy = g.Δh
        nx, ny = g.n
        V₀ = Δx * Δy
        V  = zeros(nx, ny) .+ V₀
        V[ 1,:] .*= 0.5
        V[nx,:] .*= 0.5
        V[:, 1] .*= 0.5
        V[:,ny] .*= 0.5
        return V
    end

    function cell_volume(g::AxialGrid{2})
        Δr, Δz = g.Δh
        nr, nz = g.n
        V  = zeros(nr, nz)
        for i=1:nr
            V[i,:] .= π * Δz * ((i*Δr)^2 - (i*Δr - Δr)^2)
        end
        V[ 1,:] .= π * Δz * Δr^2 * (0.75)
        V[nr,:] .= π * Δz * Δr^2 * (nr-0.25) 
        V[:, 1] .*= 0.5
        V[:,nz] .*= 0.5
        return V
    end

    function create_uniform_grid(xx, yy)
        nx, ny = length(xx), length(yy)
        xs, ys = xx[1], yy[1]
        Δx = length(xx) > 1 ? xx[2] - xx[1] : 1.0
        Δy = length(yy) > 1 ? yy[2] - yy[1] : 1.0
        x = repeat(xx,   1, ny)
        y = repeat(yy', nx,  1)
        n = nx*ny
        node = reshape(1:n, nx, ny)
        data = Dict{String,AbstractArray}()
        CartesianGrid{2}(data, (nx, ny), (Δx, Δy), node, (x, y), (xs, ys))
    end

    function create_uniform_grid(xx)
        nx = length(xx)
        xs = xx[1]
        Δx = length(xx) > 1 ? xx[2] - xx[1] : 1.0
        x = xx
        n = nx
        node = 1:n
        data = Dict{String,AbstractArray}()
        CartesianGrid{1}(data, (nx,), (Δx,), node, (x,), (xs,))
    end

    function create_axial_grid(rr, zz)
        nr, nz = length(rr), length(zz)
        rs, zs = rr[1], zz[1]
        Δr = length(rr) > 1 ? rr[2] - rr[1] : 1.0
        Δz = length(zz) > 1 ? zz[2] - zz[1] : 1.0
        r = repeat(rr,   1, nz)
        z = repeat(zz', nr,  1)
        n = nr*nz
        node = reshape(1:n, nr, nz)
        data = Dict{String,AbstractArray}()
        AxialGrid{2}(data, (nr, nz), (Δr, Δz), node, (r, z), (rs, zs))
    end

    function create_staggered_grid(g::CartesianGrid{2})
        xx, yy = g.coords
        x0, xn = extrema(xx)
        y0, yn = extrema(yy)
        Δx, Δy = g.Δh
        nx, ny = g.n
        xs = range(x0-Δx/2.0, xn+Δx/2.0, length=nx+1)
        ys = range(y0-Δy/2.0, yn+Δy/2.0, length=ny+1)
        create_uniform_grid(xs, ys)
    end

    function create_staggered_grid(g::CartesianGrid{1})
        xx, = g.coords
        x0, xn = extrema(xx)
        Δx, = g.Δh
        nx, = g.n
        xs = range(x0-Δx/2.0, xn+Δx/2.0, length=nx+1)
        create_uniform_grid(xs)
    end
end
