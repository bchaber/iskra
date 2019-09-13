using Interpolations

struct CrossSection
	nodes :: Array{Float64,2}
	interpolation 
end

CrossSection(xs :: AbstractRange, ys :: Array{Float64,1}) = CrossSection(zeros(1,1), CubicSplineInterpolation(xs, ys))
(σ :: CrossSection)(x :: Float64) = σ.interpolation(x)

Base.show(io :: IO, σ :: CrossSection) = print("σ")