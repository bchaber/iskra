    struct UniformGrid{C,D}
         n::NTuple{D,Int64}
        Δh::NTuple{D,Float64}
       bcs::NTuple{D,Tuple{Symbol, Symbol}}
     cells::Array{Int64,D}
    facets::NTuple{D,Array{Float64,D}}
    origin::NTuple{D,Float64}
    end

    CartesianGrid{D} = UniformGrid{:xy, D}

    Base.size(g::UniformGrid) = g.n
    Base.length(g::UniformGrid) = prod(g.n)

    function cell_volume(g::CartesianGrid{1})
        Δx, = g.Δh
        return Δx
    end

    function create_uniform_grid(xf; left=:periodic, right=:periodic)
        xs = first(xf)
        nx = length(xf)
        Δx = length(xf) > 1 ? xf[2] - xf[1] : 1.0
        bcs = (left, right)
        cells = collect(1:nx)
        CartesianGrid{1}(tuple(nx), tuple(Δx), tuple(bcs), cells, tuple(xf), tuple(xs))
    end