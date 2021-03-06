module ParticleInCell
  export init, solve
  export KineticSpecies, FluidSpecies
  export MaxwellianSource, DensitySource
  using Diagnostics
  using LinearAlgebra
  using Circuit
  using FiniteDifferenceMethod
  using RegularGrids

  include("pic/kinetic.jl")
  include("pic/fluid.jl")
  include("pic/cloud_in_cell.jl")
  include("pic/pushers.jl")
  include("pic/sources.jl")
  include("pic/surfaces/build.jl")
  include("pic/surfaces/track.jl")
  include("pic/surfaces/hit.jl")
  include("pic/surfaces/check.jl")
  include("pic/surfaces/wrap.jl")
  include("pic/circuit_coupling.jl")
  
  function particle_velocity(pv, p)
    vx, vy = view(pv, p, :)
    return vx, vy 
  end

  function particle_cell(px, p, Δh)
    fij = @views 1 .+ px[p, :] ./ Δh
    ij = floor.(Int64, fij)
    hxy = fij .- ij
    i, j = ij
    hx, hy = hxy
    return i, j, hx, hy
  end

  # hooks
  function enter_loop() end
  function after_loop(it, t, Δt) end
  function exit_loop() end
  function after_push(part, grid) wrap!(part, grid) end

  # extension points
  function perform!(interaction, E, Δt, config) end
  function advance!(species, E, Δt, config) end

  function init(src, species, Δt)
    sample!(src, species, Δt)
  end

  function advance!(part::KineticSpecies, E, B, Δt, config)
    tracker = config.tracker
    pusher = config.pusher
    grid = config.grid

    track!(tracker, part, Δt)
    partE = grid_to_particle(grid, part, (i,j) -> E[i,j,:])
    partB = grid_to_particle(grid, part, (i,j) -> B[i,j,:])
    push_particles!(pusher, part, partE, partB, Δt)
    check!(tracker, part, Δt)
    after_push(part, grid)

    @particle part.name*"/id" "1"  part.id part
    @particle part.name*"/mass" "kg" [part.m] part
    @particle part.name*"/charge" "C" [part.q] part
    @particle part.name*"/weighting" "1"  part.wg part weighted=true
    @particle part.name*"/momentum" "kg*m/s" part.m*part.v  part withcomponents=true
    @particle part.name*"/position" "m"   part.x  part withcomponents=true
    @particle part.name*"/positionOffset/x" "m" [0.] part
    @particle part.name*"/positionOffset/y" "m" [0.] part
    @particle part.name*"/positionOffset/z" "m" [0.] part
  end
  
  function advance!(fluid :: FluidSpecies, E, B, Δt, config)
    q, m = fluid.q, fluid.m
    v = E*Δt*(q/m)
    Δn = calculate_advection_diffusion(fluid.n, fluid.μ, v, config.grid.Δh, Δt)
    fluid.n .+= Δn
    
    @field "d"*fluid.name "1/m^2" Δn config.grid
    @field "v"*fluid.name "m/s"  v config.grid withcomponents=true
  end

  function solve(config, Δt=1e-5, timesteps=200)
    pusher = config.pusher
    interactions = config.interactions
    sources = config.sources
    species = config.species
    circuit = config.circuit
    solver = config.solver
    cells = config.cells
    grid  = config.grid
    Δt = float(Δt)

    enter_loop()

    ϕ = zeros(size(grid))
    ρ = zeros(size(grid))
    E = zeros(size(grid)..., 3)
    B = zeros(size(grid)..., 3)

    for iteration in 1:timesteps # iterate for ts time step
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
        advance!(part, E, B, Δt, config)
      end
      advance!(circuit, ϕ, Δt, config)
      # Calculate charge density
      fill!(ρ, 0.0)
      for part in species
        part.n = density(part, grid)
        ρ .+= part.n .* part.q

        @field "n"*part.name "1/m^2" part.n grid
      end
      # Calculate electric field
      ϕ = calculate_electric_potential(solver, -ρ)
      E = calculate_electric_field(solver, ϕ)
      B = calculate_magnetic_field(solver)

      @field "rho" "C/m^2" ρ grid
      @field "phi" "V"     ϕ grid
      @field "E"   "V/m"   E grid withcomponents=true
      @field "B"   "A/m"   B grid withcomponents=true
      after_loop(iteration, iteration*Δt-Δt, Δt)
    end

    exit_loop()
    println("Complete!")
  end
end
