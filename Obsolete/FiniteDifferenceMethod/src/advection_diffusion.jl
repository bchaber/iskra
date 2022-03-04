struct AdvectionDiffusionSolver
end

function calculate_advection_diffusion(n::Vector{Float64}, D, v, Δh, Δt)
    nx, = size(n)
    Δx, = Δh
    vx = view(v,:,1)
    Δn = zeros(nx)

    for i=2:nx-1
      # diffusion
      Δn[i] += Δt*D*(n[i-1] - 2n[i] + n[i+1])/Δx^2
      # assume no advection
      #Δn[i] -= Δt*n[i] * (vx[i+1] - vx[i-1])/2Δx
      #Δn[i] -= Δt*(vx[i] * (n[i+1] - n[i-1]))/Δx 
    end

    return Δn
end

function calculate_advection_diffusion(n::Matrix{Float64}, D, v, Δh, Δt)
    nx, ny = size(n)
    Δx, Δy = Δh
    vx = view(v,:,:,1)
    vy = view(v,:,:,2)
    Δn = zeros(nx, ny)

    for i=2:nx-1
      for j=2:ny-1
        # diffusion
        Δn[i,j] += Δt*D*(n[i-1,j] + n[i,j-1] - 4n[i,j] + n[i+1,j] + n[i,j+1])/(Δx*Δy)
        # advection
        Δn[i,j] -= Δt*n[i,j] * (vx[i+1,j] - vx[i-1,j] + vy[i,j+1] - vy[i,j-1])/2Δx
        Δn[i,j] -= Δt*(vx[i,j] * (n[i+1,j] - n[i-1,j]) + vy[i,j] * (n[i,j+1] - n[i,j-1]))/Δx 
      end
    end

    return Δn
end