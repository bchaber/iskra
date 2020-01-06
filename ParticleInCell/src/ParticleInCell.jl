module ParticleInCell
  using LinearAlgebra
  using Circuit
  using FiniteDifferenceMethod

  include("pic/kinetic.jl")
  include("pic/fluid.jl")
  include("pic/cloud_in_cell.jl")
  include("pic/diagnostics.jl")
  include("pic/pushers.jl")
  include("pic/sources.jl")
  include("pic/surface.jl")
  
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

  # extension points
  function perform!(interaction, E, Δt, config) end
  function advance!(species, E, Δt, config) end

  function init(src, species, Δt)
    sample!(src, species, Δt)
  end

  function advance!(circuit :: Nothing, ϕ, Δt, config, ε0) end
  function advance!(circuit :: CircuitRLC, ϕ, Δt, config, ε0)
    nx, ny = size(config.grid)
    Δx, Δy, Δz = config.grid.Δh
    A = ny * Δy*Δz
    V = ϕ[1,1] - ϕ[nx,ny]
    advance_circuit!(circuit, V, Δt)
    dσ = -Δt*circuit.i/A
    @diag "dσ" TimeData(dσ)
    config.solver.b[config.solver.dofs[:σ][1]] += dσ/ε0
  end

  function advance!(part :: KineticSpecies, E, Δt, config)
    tracker = config.tracker
    pusher = config.pusher
    grid = config.grid
    nx, ny = size(grid)

    track!(tracker, part)
    partE = grid_to_particle(grid, part, (i,j) -> E[i, j, :])
    push_particles!(pusher, part, partE, Δt)
    check!(tracker, part, Δt)

    @diag "pv"*part.name ParticleVectorData(part.x,part.v,part.id, part.wg, part.np)
    @diag "pE"*part.name ParticleVectorData(part.x,partE, part.id, part.wg, part.np)
  end
  
  function advance!(fluid :: FluidSpecies, E, Δt, config)
    q, m = fluid.q, fluid.m
    v = E*Δt*(q/m)
    Δn = calculate_advection_diffusion(fluid.n, fluid.μ, v, config.grid.Δh, Δt)
    fluid.n .+= Δn
    
    @diag  "v"*fluid.name GridData( v, config.grid.x, config.grid.y)
  end

  function solve(config, Δt=1e-5, timesteps=200, ε0=1.0)
    pusher = config.pusher
    interactions = config.interactions
    sources = config.sources
    species = config.species
    circuit = config.circuit
    solver = config.solver
    cells = config.cells
    grid  = config.grid
    
    nx, ny = size(grid)
    Δx, Δy, Δz = grid.Δh
    spacing = [Δx, Δy]
    origin  = grid.origin
    enter_loop()

    ϕ = zeros(nx, ny)
    ρ = zeros(nx, ny)
    E = zeros(nx, ny, 3)

    for iteration=1:timesteps # iterate for ts time step
      # Create particles
      for src in sources
        sample!(src, src.species, Δt)
      end

      # Solve reactions
      for interaction in interactions
        perform!(interaction, E, Δt, config)
      end
      # Advance species
      for part in species
        advance!(part, E, Δt, config)
      end
      advance!(circuit, ϕ, Δt, config, ε0)
      # Calculate charge density
      fill!(ρ, 0.0)
      for part in species
        part.n = density(part, grid)
        ρ .+= part.n .* part.q

        @diag "n"*part.name NodeData(part.n, origin, spacing)
      end
      # Calculate electric field
      ϕ  = calculate_electric_potential(solver, -ρ/ε0)
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
