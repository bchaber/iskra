module RegularGrid
    export create_uniform_grid, create_staggered_grid
    export UniformGrid, CartesianGrid, AxialGrid

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
