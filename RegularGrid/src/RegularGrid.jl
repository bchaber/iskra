module RegularGrid
    export create_uniform_grid, create_staggered_grid
    export UniformGrid, CoordinateSystem, XY2D, RZ2D
    abstract type CoordinateSystem end
    abstract type XY2D <: CoordinateSystem end
    abstract type RZ2D <: CoordinateSystem end
    struct UniformGrid{cs <: CoordinateSystem}
      data::Dict{String,AbstractArray}
        nx::Integer
        ny::Integer
        nz::Integer
        Δh::Tuple{Float64,Float64,Float64}
         x::AbstractArray
         y::AbstractArray
         z::AbstractArray
      node::AbstractArray
    origin::AbstractArray
    end

    Base.size(g::UniformGrid{<:CoordinateSystem}) = (g.nx, g.ny)
    Base.length(g::UniformGrid{<:CoordinateSystem}) = g.nx*g.ny
    Base.getindex(g::UniformGrid{<:CoordinateSystem}, k) = g.data[k]
    Base.setindex!(g::UniformGrid{<:CoordinateSystem}, v, k) =
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
        UniformGrid{XY2D}(data, nx, ny, 0, (Δx, Δy, 1.), x, y, z, node, [xs,ys,0.])
    end

    function create_staggered_grid(g::UniformGrid{XY2D})
        x0, xn = extrema(g.x)
        y0, yn = extrema(g.y)
        Δx, Δy, ~ = g.Δh
        xs = range(x0-Δx/2, xn+Δx/2, length=g.nx+1)
        ys = range(y0-Δy/2, yn+Δy/2, length=g.ny+1)
        create_uniform_grid(xs, ys)
    end
end
