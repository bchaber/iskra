module RegularGrid
    export create, UniformGrid
    struct UniformGrid
        nx::Integer
        ny::Integer
        Δh::AbstractFloat
        x::AbstractArray
        y::AbstractArray
        dof::AbstractArray
        origin::AbstractArray
    end

    Base.size(g::UniformGrid) = (g.nx, g.ny)
    Base.length(g::UniformGrid) = g.nx*g.ny

    function UniformGrid(xs, ys, xe, ye, Δh)
        xx, yy = xs:Δh:xe, ys:Δh:ye
        UniformGrid(xx, yy)
    end
   
    function UniformGrid(xx, yy)
        nx, ny = length(xx), length(yy)
        xs, ys = xx[1], yy[1]
        Δh = xx[2] - xx[1]
        x = repeat(xx, 1, ny)
        y = repeat(yy, 1, nx)'
        n = nx*ny
        dof = reshape(1:n, nx, ny)
        return UniformGrid(nx, ny, Δh, x, y, dof, [xs,ys])
    end
end
