module ParticleInCell
  using FiniteDifferenceMethod
  
  include("sugar.jl")

  include("pic/kinetic.jl")
  include("pic/fluid.jl")
  include("pic/cloud_in_cell.jl")
  include("pic/diagnostics.jl")
  include("pic/pushers.jl")
  include("pic/sources.jl")

  function particle_cell(px, p, Δh)
    fij = @views 1 .+ px[p, :] ./ Δh
    ij = floor.(Int64, fij)
    hxy = fij .- ij
    i, j  = ij
    hx,hy = hxy
    return i, j, hx, hy
  end

  # hooks
  function enter_loop() end
  function after_loop(it) end
  function exit_loop() end

  function init(src, Δt)
    sample!(src, src.species, Δt)
  end

  function advance!(part :: KineticSpecies, E, Δt, config)
    pusher = config.pusher
    grid = config.grid
    nx, ny = size(grid)

    partE = grid_to_particle(grid, part, (i,j) -> E[i, j, :])
    push_particles!(pusher, part, partE, Δt)
    remove_particles!(part, grid.Δh, (i,j) -> i < 1 || i >= nx || j < 1 || j >= ny)

    @diag "pv"*part.name ParticleVectorData(part.x,part.v,part.id, part.np)
    @diag "pE"*part.name ParticleVectorData(part.x,partE, part.id, part.np)
  end
  
  function advance!(fluid :: FluidSpecies, E, Δt, config)
    q, m = fluid.q, fluid.m
    v = E*Δt*(q/m)
    Δn = calculate_advection_diffusion(fluid.n, v, config.grid.Δh, Δt)
    fluid.n .+= Δn  

    @diag "Δn"*fluid.name NodeData(Δn, config.grid.origin, [1,1]*config.grid.Δh)
    @diag  "v"*fluid.name GridData( v, config.grid.x, config.grid.y)
  end

  function solve(config, Δt=1e-5, timesteps=200, ε0=1.0)
    pusher = config.pusher
    sources = config.sources
    species = config.species
    solver = config.solver
    cells = config.cells
    grid  = config.grid
    
    nx, ny = size(grid)
    Δh = grid.Δh
    spacing = [1,1]*grid.Δh
    origin  = grid.origin
    enter_loop()

    ρ = zeros(nx, ny)
    E = zeros(nx, ny, 3)

    @diag "εr" GridData(cells["εr"], cells.x,cells.y)

    for iteration=1:timesteps # iterate for ts time step
      # Create particles
      for src in sources
        sample!(src, src.species, Δt)
      end
      # Advance species
      for part in species
        advance!(part, E, Δt, config)
      end
      # Calculate charge density
      fill!(ρ, 0.0)
      for part in species
        n = density(part, grid)
        ρ .+= n .* part.q

        @diag "n"*part.name NodeData(n, origin, spacing)
      end
      # Calculate electric field
      Q  = ρ .* Δh .^2
      ϕ  = calculate_electric_potential(solver, -Q/ε0)
      E  = calculate_electric_field(solver, ϕ)

      @diag "ρ" NodeData(ρ, origin, spacing)
      @diag "ϕ" NodeData(ϕ, origin, spacing)
      @diag "E" GridData(E, grid.x,  grid.y)

      after_loop(iteration)

      println("Time Step #", iteration)#, ", Particles #", [part.np for part in config.species])
    end

    exit_loop()
    println("Complete!")
  end
end
