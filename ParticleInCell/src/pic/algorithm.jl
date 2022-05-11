
  function particle_cell(part::KineticSpecies{1, V}, p, Δh) where {V}
    fij = 1.0 .+ part.x[p] ./ Δh
    ij = floor.(Int64, fij)
    hxy = fij .- ij
    i, = ij
    hx, = hxy
    return i, 1, hx, 0.
  end

  function number_density!(n, species :: KineticSpecies, grid)
    n .= 0.0
    particle_to_grid!(n, species, grid, species.wg)
    n ./= cell_volume(grid)
    return nothing
  end

  function advance!(fluid::FluidSpecies, E, B, Δt, pusher, grid) end
  function advance!(part::KineticSpecies, E, B, Δt, pusher, grid)
    data = pusher.data[part.name]
    grid_to_particle!(data.E, grid, part, E)
    grid_to_particle!(data.B, grid, part, B)
    push_particles!(pusher, part, Δt)
    after_push(part, grid)
    return nothing
  end
  
  function average_electric_field(E)
    for i = 2:length(E)
      E[i-1] += E[i]
      E[i-1] /= 2.0
    end
    return nothing
  end

  function current_deposition(j, part, Δx, Δt)
    j .= 0.0
    c  = Δx/Δt
    nx = length(j) - 1
    for p=1:part.np
      i, _, h, _ = particle_cell(part, p, grid.Δh)
      vp = part.v[p][1]
      hh = h - vp / c
      di = fld(hh, 1.0)
      
      h1 = mod(hh, 1.0)
      h2 = h

      b1 = abs(0.5 - h1)
      b2 = abs(0.5 - h2)
      
      if vp > 0.0
        #print("vp > 0")
        if di == 0
          #print(" i1=i2")
          if     h1 < 0.5 && h2 > 0.5 #println(" a)")
            j[i]   += c * b1
            j[i+1] += c * b2
            if i == nx
              j[1] += c * b2
            end
          elseif h1 > 0.5 && h2 > 0.5 #println(" b)")
            j[i+1]   += c * -(b1 - b2)
            if i == nx
              j[1] += c * -(b1 - b2)
            end
          elseif h1 < 0.5 && h2 < 0.5 #println(" c)")
            j[i] += c * -(b2 - b1)
            if i == 1
              j[nx + 1] += c * -(b2 - b1)
            end
          end

        else # moves to another cell
          #print(" i1/i2")
          if     h1 > 0.5 && h2 < 0.5 #println(" a)")
            j[i]   += c * (1. - b1 - b2)
            if i == 1
              j[nx + 1] += c * (1. - b1 - b2)
            end
          elseif h1 > 0.5 && h2 > 0.5 #println(" b)")
            j[i]   += c * (1. - b1)
            j[i+1] += c * b2
          elseif h1 < 0.5 && h2 < 0.5 #println(" c)")
            j[i > 1 ? i-1 : nx + 1] += c * b1
            j[i]   += c * (1. - b2)
          end

        end

      else # vp < 0.0
        #print("vp < 0")
        if di == 0
          #print(" i1=i2")
          if     h1 > 0.5 && h2 < 0.5 #println(" a)")
            j[i]   += c * (-b2)
            j[i+1] += c * (-b1)
          elseif h1 > 0.5 && h2 > 0.5 #println(" b)")
            j[i+1]   += c * -(b1 - b2)
            if i == nx
              j[1] += c * -(b1 - b2)
            end
          elseif h1 < 0.5 && h2 < 0.5 #println(" c)")
            j[i] += c * -(b2 - b1)
            #if i == 1
            #  j[nx + 1] += c * -(b2 - b1)
            #end
          end

        else # moves to another cell
          #print(" i1/i2")
          if     h1 < 0.5 && h2 > 0.5 #println(" a)")
            j[i+1] += c * -(1. - b1 - b2)
            if i == nx
              j[1] += c * -(1. - b1 - b2)
            end
          elseif h1 > 0.5 && h2 > 0.5 #println(" b)")
            j[i < nx ? i+2 : 2] += c *-(1. - b2)
            j[i+1] += c * (+b1)
          elseif h1 < 0.5 && h2 < 0.5 #println(" c)")
            j[i]   += c * (-b2)
            j[i+1] += c * -(1. - b1)
          end

        end

      end
    end
    return nothing
  end

function solve(config, Δt=1e-5, timesteps=200)
    diagnostics = config.diagnostics
    pusher = config.pusher
    interactions = config.interactions
    species = config.species
    solver = config.solver
    grid  = config.grid
    Δt = float(Δt)

    ϕ = solution(solver)
    ∇ = gradient(solver)
    n = zeros(Float64, size(grid) .- 1)
    j = zeros(Float64, size(grid))
    ρ = zeros(Float64, size(grid) .- 1)
    E = zeros(SVector{3, Float64}, size(grid)...)
    B = zeros(SVector{3, Float64}, size(grid)...)
    Jx = zeros(length(grid))
    Ex = view(reinterpret(Float64, E), 1:3:3length(E))
    Ey = view(reinterpret(Float64, E), 2:3:3length(E))
    Ez = view(reinterpret(Float64, E), 3:3:3length(E))
    enter_loop()
    for iteration in 1:timesteps # iterate for ts time step
      # Calculate charge number density
      fill!(ρ, 0.0)
      for part in species
        number_density!(n, part, grid)
        ρ .+= n * part.q
        #diagnostics["n"*part.name][:, iteration] .= n
      end
      # Calculate electric field
      calculate_electric_potential(solver, copy(ρ))
      ∇(ϕ; result=E)
      average_electric_field(E)
      #if 700 > iteration > 100
      #  v = @SVector [1.5e9sin(2pi*iteration/100.0), 0.0, 0.0]
      #  for j=30:33 E[j] = v end
      #  for j=90:93 E[j] = v end
      #end
      # Solve reactions
      for interaction in interactions
        perform!(interaction, E, Δt, config)
      end
      # Advance species
      for part in species
        advance!(part, E, B, Δt, pusher, grid)
      end
      # Calculate current density
      fill!(Jx, 0.0)
      #print("it. ", iteration, " ")
      for part in species
        current_deposition(j, part, first(grid.Δh), Δt)
        #println(part, " ", maximum(j))
        Jx .+= j * part.q
      end
      diagnostics["rho"][:, iteration] .= ρ
      diagnostics["phi"][:, iteration] .= ϕ
      diagnostics["Ex"][:,  iteration] .= Ex
      diagnostics["Jx"][:,  iteration] .= Jx
      after_loop(iteration, iteration*Δt-Δt, Δt)
    end
    exit_loop()
    println("Complete!")
  end
