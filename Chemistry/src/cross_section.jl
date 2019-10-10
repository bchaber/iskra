using Interpolations

struct CrossSection
	nodes :: Array{Float64,1}
	interpolation
end

CrossSection(xs :: AbstractRange, ys :: Array{Float64,1}) =
	CrossSection(ys, CubicSplineInterpolation(xs, ys; extrapolation_bc=.0))
(σ :: CrossSection)(x :: Float64) = max(0, σ.interpolation(x))
Base.maximum(σ :: CrossSection) = maximum(σ.nodes)
Base.show(io :: IO, σ :: CrossSection) = print(io, "σ: ")