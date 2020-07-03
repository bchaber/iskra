module RegularGrid
    export create_uniform_grid, create_staggered_grid
    export UniformGrid, CartesianGrid, AxialGrid

    struct UniformGrid{C, D}
      data::Dict{String,AbstractArray}
         n::NTuple{D,Integer}
        Δh::NTuple{D,Float64}
         u::Array{Float64,D}
         v::Array{Float64,D}
         w::Array{Float64,D}
      node::Array{Float64,D}
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
        z = zeros(nx, ny)
        n = nx*ny
        node = reshape(1:n, nx, ny)
        data = Dict{String,AbstractArray}()
        CartesianGrid{2}(data, (nx, ny), (Δx, Δy), x, y, z, node, (xs, ys))
    end

    function create_staggered_grid(g::UniformGrid{:xy, 2})
        x0, xn = extrema(g.u)
        y0, yn = extrema(g.v)
        Δx, Δy = g.Δh
        nx, ny = g.n
        xs = range(x0-Δx/2, xn+Δx/2, length=nx+1)
        ys = range(y0-Δy/2, yn+Δy/2, length=ny+1)
        create_uniform_grid(xs, ys)
    end
end
