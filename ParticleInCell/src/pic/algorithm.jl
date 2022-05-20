
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

  function after_solver(state, config, iteration) end
  function after_pusher(state, config, iteration) end
  
  function advance!(fluid::FluidSpecies, E, B, Δt, pusher, grid) end
  function advance!(part::KineticSpecies, E, B, Δt, pusher, grid)
    data = pusher.data[part.name]
    grid_to_particle!(data.E, grid, part, E)
    grid_to_particle!(data.B, grid, part, B)
    push_particles!(pusher, part, Δt)
    return nothing
  end
  
  function average_over_cells(F)
    for i = 2:length(F)
      F[i-1] += F[i]
      F[i-1] /= 2.0
    end
    return nothing
  end

  function current_deposition!(j, part, grid, c)
    j .= 0.0
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
            #if i == nx
            #  j[1] += c * -(b1 - b2)
            #end
          elseif h1 < 0.5 && h2 < 0.5 #println(" c)")
            j[i] += c * -(b2 - b1)
            #if i == 1
            #  j[nx + 1] += c * -(b2 - b1)
            #end
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
            if i > 1
              j[i-1]    += c * b1
            else
              j[nx + 1] += c * b1
            end
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
            if i < nx
              j[i+2] += c *-(1. - b2)
            else
              j[2]   += c *-(1. - b2)
            end
            j[i+1] += c * (+b1)
          elseif h1 < 0.5 && h2 < 0.5 #println(" c)")
            j[i]   += c * (-b2)
            j[i+1] += c * -(1. - b1)
          end

        end

      end
    end
    j .*= part.w0
    j ./= cell_volume(grid)
    return nothing
  end

function solve(state, pusher, solver, grid, interactions, timesteps, context)
    species = state.particles
    Δt = state.timestep
    Δx = state.cellsize

    c = Δx / Δt
    ϕ = solution(solver)
    ∇ = gradient(solver)

    ρ = state.density
    E = state.electric
    B = state.magnetic
    J = state.current
    
    n = similar(ρ)
    j = similar(J)
    
    Ex = view(reinterpret(Float64, E), 1:3:3length(E))
    Ey = view(reinterpret(Float64, E), 2:3:3length(E))
    Ez = view(reinterpret(Float64, E), 3:3:3length(E))
    Bx = view(reinterpret(Float64, B), 1:3:3length(B))
    By = view(reinterpret(Float64, B), 2:3:3length(B))
    Bz = view(reinterpret(Float64, B), 3:3:3length(B))
    
    enter_loop(context, state)
    for iteration in 1:timesteps # iterate for ts time step
      t = Δt * (iteration - 1)
      # Calculate charge number density
      fill!(ρ, 0.0)
      for part in species
        number_density!(n, part, grid)
        ρ .+= n * part.q
      end
      # Calculate electric field
      calculate_electric_potential(solver, copy(ρ))
      ∇(ϕ; result=E)
      average_over_cells(E)
      average_over_cells(J)
      after_solver(context, state, iteration)
      
      # Solve reactions
      for interaction in interactions
        perform!(interaction, E, Δt, config)
      end
      # Advance species
      for part in species
        advance!(part, E, B, Δt, pusher, grid)
      end
      # Calculate current density
      after_pusher(context, state, iteration)
      fill!(J, 0.0)
      for part in species
        current_deposition!(j, part, grid, c)
        J .+= j * part.q
      end

      after_loop(context, state, iteration)
    end
    exit_loop(context, state)
    println("Complete!")
  end
  
