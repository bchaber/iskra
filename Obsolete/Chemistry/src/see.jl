# Going to be implemented according to 
# D. Sydorenko, "Particle-in-cell simulations of electron dynamics in low pressure
# discharges with magnetic fields,” 2006, PhD. Thesis, University of Saskatchewan, Canada,
# url: http://hdl.handle.net/10388/etd-06142006-111353 
# and extended to ion-induced secondary electron emission

using ParticleInCell
using Distributions
using LinearAlgebra
function true_secondary_energy(part::KineticSpecies{D,V}) where {D,V}
	μ, σ = 1.65, 1.1
	rand(LogNormal(μ, σ))
end

# Vaughan coefficient
function vaughan(;w₀, w₀max, γ₀max, kₛ=0.0)
	wmax(θ) = w₀max*(1.0 + kₛ/1π * θ^2)
	γmax(θ) = γ₀max*(1.0 + kₛ/2π * θ^2)
	v(w, θ) = w > w₀ ? (w - w₀)/(wmax(θ) - w₀) : 0.0

	function (w, θ)
		k = w > wmax(θ) ? 0.25 : 0.62
		γmax(θ)*(v(w,θ)*exp(1.0 - v(w,θ)))^k
	end
end

# elastic reflection coefficient
function elastic(γᵥ::Function, wₑ, wₑmax, γₑmax; Δₑ=13., rₑ=0.03)
	v₁(w) = (w - wₑ) / (wₑmax - wₑ)
	v₂(w) = (w - wₑmax) / Δₑ

	function (w, θ)
		if wₑ < w <= wₑmax
			rₑ*γᵥ(w, θ) + γₑmax*v₁(w)*exp(1.0 - v₁(w))
		elseif   w > wₑmax
			rₑ*γᵥ(w, θ) + γₑmax*(1.0 + v₂(w))*exp(-v₂(w))
		else
			0.0
		end
	end
end

# inelastic reflection coefficient
function inelastic(γᵥ::Function; rᵢ=0.07)
	function (w, θ)
		rᵢ*γᵥ(w, θ)
	end
end

# true secondary emission
function secondary(γᵥ::Function; rₑ, rᵢ)
	function (w, θ)
		(1.0 - rₑ - rᵢ)*γᵥ(w, θ)
	end
end

γᵥ = vaughan(w₀=13., w₀max=500., γ₀max=3., kₛ=1.)
γₑ = elastic(γᵥ, 2., 10., 0.55, rₑ=0.03)
γᵢ = inelastic(γᵥ, rᵢ=0.07)
γₜ = secondary(γᵥ, rₑ=0.03, rᵢ=0.07)

@inline mass(species) = species.m/1.60217646e-19 # qe = 1.60217646e-19C

@inline both(α)  = α ≠ 0.0
@inline lower(α) = α < 0.0
@inline upper(α) = α > 0.0

function predicate(boundary)
	if boundary == :left   return lower end
	if boundary == :right  return upper end
	if boundary == :top    return upper end
	if boundary == :bottom return lower end
	if boundary == :back   return lower end
	if boundary == :front  return upper end
	return both
end

function normal(boundary)
	if boundary == :left   return [-1., 0., 0.] end
	if boundary == :right  return [ 1., 0., 0.] end
	if boundary == :top    return [ 0., 1., 0.] end
	if boundary == :bottom return [ 0.,-1., 0.] end
	if boundary == :back   return [ 0., 0.,-1.] end
	if boundary == :front  return [ 0., 0., 1.] end
	return [0.,0.,0]
end

function dims(boundary, D)
	if boundary == :left   return (1,) end
	if boundary == :right  return (1,) end
	if boundary == :top    return (2,) end
	if boundary == :bottom return (2,) end
	if boundary == :back   return (3,) end
	if boundary == :front  return (3,) end
	return 1:D
end

γ₀(w, θ) = 0.0

snells_law(vᵢ, n̂) = vᵢ .- 2.0(n̂ ⋅ vᵢ)*n̂	

function inject_secondary!(secondary::KineticSpecies{D,V}, x, n̂, dt) where {D,V}
	m = mass(secondary)
	np = secondary.np
	px = view(secondary.x, np+1,:)
	pv = view(secondary.v, np+1,:)
	ε = true_secondary_energy(secondary)
	v = diffuse_reflection(n̂) * sqrt(2.0ε/m)
	pv[1:V] .= v
	px[1:D] .= dt * v[1:D] .+ x[1:D]
	secondary.np += 1
end

function emit!(primary::KineticSpecies{D,V}, secondary::KineticSpecies{D,V}, 
	grid, material::Symbol; boundary=:all,
	 γₜ, γₑ=γ₀, γᵢ=γ₀) where {CS, D, V}
	hit = predicate(boundary)
	n̂ = normal(boundary)

	for i in dims(boundary, D)
	  nx = grid.n[i] - 1
	  Δx = grid.Δh[i]
	  Lx = nx*Δx
	  ox = grid.origin[i]
	  px = view(primary.x, 1:primary.np, i)
	  for p=reverse(1:primary.np)
	    α = fld(px[p] - ox, Lx)
	    if hit(α)
	    	pv = view(primary.v, p, 1:V)
	    	vₚ = view(primary.v, p, 1:D)
	    	mₚ = mass(primary)
	    	θ = sinθ = norm(pv × n̂)/norm(pv)
	    	w = 0.5mₚ * dot(pv,pv)
	    	R₁ = rand()
	    	γe = γₑ(w, θ)
	    	γi = γᵢ(w, θ)
	    	γt = γₜ(w, θ)

	    	xₚ = view(primary.x, p, :)
	    	dt = mod(px[p], Lx) / abs(pv[i])
	    	#println("hit! ", p, " ", xₚ, " (", xₚ .- vₚ .* dt, ") ", w, "eV")
	    	#println(γe, "/", γi, "/", γt, ":", R₁)

	    	x0 = xₚ .- vₚ * dt
	    	if γe + γi > R₁ > γe
	    		xₚ .-= vₚ * dt
	    		pv  .= rand()*snells_law(pv, n̂) # specular
	    		#pv .= rand()*norm(pv)*diffuse_reflection(n̂) # diffuse
	    		xₚ .+= vₚ * dt
	    		#println("inelastic collision")
	    		continue
	    	end

	    	if R₁ < γe
	    		xₚ .-= vₚ * dt
	    		pv  .= snells_law(pv, n̂) # specular
	    		#pv .= norm(pv)*diffuse_reflection(n̂) # diffuse
	    		xₚ .+= vₚ * dt
	    		#println("elastic collision")
	    		continue
	    	end
	    	
	    	γ = γe + γi + γt
	    	while γ > 1.0
	    		inject_secondary!(secondary, x0, n̂, dt) # inject true secondary
	    		println("secondary!")
	    		γ = γ - 1.0
	    	end

	    	R₂ = rand()
	    	if R₂ < γ
	    		inject_secondary!(secondary, x0, n̂, dt) # inject true secondary
	    		println("secondary!!!")
	    	else
	    		ParticleInCell.remove!(primary, p) # absorb
	    		#println("absorbed")
	    	end
	    end
	  end
	end
end