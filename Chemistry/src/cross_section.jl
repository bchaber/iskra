using Interpolations

struct CrossSection
	nodes :: Array{Float64,2}
	interpolation 
end

CrossSection(xs :: AbstractRange, ys :: Array{Float64,1}) =
	CrossSection(zeros(1,1), CubicSplineInterpolation(xs, ys; extrapolation_bc=.0))
(σ :: CrossSection)(x :: Float64) = max(0, σ.interpolation(x))

Base.show(io :: IO, σ :: CrossSection) = print(io, "σ: ")