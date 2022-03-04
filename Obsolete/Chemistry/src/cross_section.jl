using Interpolations

struct CrossSection
	nodes :: Array{Float64,2}
	interpolation
end

CrossSection(xs :: AbstractVector{Float64}, ys :: AbstractVector{Float64}) =
	CrossSection([xs ys], LinearInterpolation(xs, ys;
		extrapolation_bc=Flat()))
CrossSection(nodes :: Array{Float64,2}) =
	CrossSection(nodes, LinearInterpolation(nodes[:,1], nodes[:,2];
		extrapolation_bc=Flat()))
(σ :: CrossSection)(x :: Float64) = σ.interpolation(x)		
Base.maximum(σ :: CrossSection) = maximum(σ.nodes[:,2])
Base.argmax(σ :: CrossSection) = σ.nodes[argmax(σ.nodes[:,2]),1]
Base.show(io :: IO, σ :: CrossSection) = print(io, "σ: ")