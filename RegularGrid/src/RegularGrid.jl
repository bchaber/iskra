module RegularGrid
    export create, UniformGrid
    struct UniformGrid
      data::Dict{String,AbstractArray}
        nx::Integer
        ny::Integer
        Δh::Tuple{Float64,Float64,Float64}
         x::AbstractArray
         y::AbstractArray
         z::AbstractArray
      node::AbstractArray
    origin::AbstractArray
    end

    Base.size(g::UniformGrid) = (g.nx, g.ny)
    Base.length(g::UniformGrid) = g.nx*g.ny
    Base.getindex(g::UniformGrid, k) = g.data[k]
    Base.setindex!(g::UniformGrid, v, k) = if length(v) ≠ length(g) error("Size mismatch!") else setindex!(g.data, v, k) end

    function create_uniform_grid(xx, yy)
        nx, ny = length(xx), length(yy)
        xs, ys = xx[1], yy[1]
        Δx = xx[2] - xx[1]
        Δy = yy[2] - yy[1]
        Δz = 1.0
        x = repeat(xx,   1, ny)
        y = repeat(yy', nx,  1)
        z = zeros(nx, ny)
        n = nx*ny
        node = reshape(1:n, nx, ny)
        data = Dict{String,AbstractArray}()
        UniformGrid(data, nx, ny, (Δx, Δy, Δz), x, y, z, node, [xs,ys])
    end

    function create_staggered_grid(g::UniformGrid)
        x0, xn = extrema(g.x)
        y0, yn = extrema(g.y)
        Δx, Δy, ~ = g.Δh
        xs = range(x0-Δx/2, xn+Δx/2, length=g.nx+1)
        ys = range(y0-Δy/2, yn+Δy/2, length=g.ny+1)
        create_uniform_grid(xs, ys)
    end
end
