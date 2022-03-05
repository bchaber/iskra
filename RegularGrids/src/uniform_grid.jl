    struct UniformGrid{C, D}
      data::Dict{String,AbstractArray}
         n::NTuple{D,Integer}
        Δh::NTuple{D,Float64}
       bcs::NTuple{D,Tuple{Symbol, Symbol}}
      node::Array{Float64,D}
    coords::NTuple{D,Array{Float64,D}}
    origin::NTuple{D,Float64}
    end

    CartesianGrid{D} = UniformGrid{:xy, D}

    Base.size(g::UniformGrid) = g.n
    Base.length(g::UniformGrid) = prod(g.n)
    Base.getindex(g::UniformGrid, k) = g.data[k]
    Base.setindex!(g::UniformGrid, v, k) =
        length(v) ≠ length(g) ? error("Size mismatch!") : setindex!(g.data, v, k)

    function cell_volume(g::CartesianGrid{1})
        Δx, = g.Δh
        return Δx
    end

    function create_uniform_grid(xx; left=:periodic, right=:periodic)
        nx = length(xx)
        xs = xx[1]
        Δx = length(xx) > 1 ? xx[2] - xx[1] : 1.0
        x = xx
        n = nx
        bcs = (left, right)
        node = collect(1:n)
        data = Dict{String,AbstractArray}()
        CartesianGrid{1}(data, (nx,), (Δx,), (bcs,), node, (x,), (xs,))
    end