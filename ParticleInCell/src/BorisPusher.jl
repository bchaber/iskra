module BorisPusher

  export push!
  function push!(part, E, Δh, Δt)
    np = part.np
    F = part.q*E[1:np,:] # Lorentz force, F=qE
    a = F./part.m        # acceleration
    for p=1:part.np
      part.v[p,:] .= part.v[p,:] .+ Δt*a[p,:]
      part.x[p,:] .= part.x[p,:] .+ Δt*part.v[p,:]
    end
  end
end
