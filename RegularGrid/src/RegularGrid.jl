module RegularGrid
    export create, UniformGrid
    struct UniformGrid
        nx::Integer
        ny::Integer
        Δh::AbstractFloat
         x::AbstractArray
         y::AbstractArray
         z::AbstractArray
       dof::AbstractArray
    origin::AbstractArray
    end

    Base.size(g::UniformGrid) = (g.nx, g.ny)
    Base.length(g::UniformGrid) = g.nx*g.ny

    function create_uniform_grid(xx, yy)
        nx, ny = length(xx), length(yy)
        xs, ys = xx[1], yy[1]
        Δh = xx[2] - xx[1]
        x = repeat(xx,   1, ny)
        y = repeat(yy', nx,  1)
        z = zeros(nx, ny)
        n = nx*ny
        dof = reshape(1:n, nx, ny)
        return UniformGrid(nx, ny, Δh, x, y, z, dof, [xs,ys])
    end
end
