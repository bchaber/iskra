struct AdvectionDiffusionSolver
end

function calculate_advection_diffusion(n, D, v, Δh, Δt)
    nx, ny = size(n)
    Δx, Δy, ~ = Δh
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