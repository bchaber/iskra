module ParticleInCell
  export solve

  import Printf:@sprintf,println

  include("sugar.jl")

  include("pic/macroparticle.jl")
  include("pic/cloud_in_cell.jl")

  include("Pusher.jl")
  include("Source.jl")
  using .Pusher
  using .Source
  using FiniteDifferenceMethod

  function particle_cell(px, p, Δh)
    fij = 1 .+ px[p,:] ./ Δh
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

  function solve(config, Δt=1e-5, timesteps=200)
    pusher = config.pusher
    source = config.source
    species = config.species
    solver = config.solver
    grid = config.grid

    nx, ny = size(grid)
    Δh = grid.Δh
    Δh² = Δh^2
    spacing = [1,1]*grid.Δh
    origin  = grid.origin
    part = config.species[1]

    enter_loop()

    for iteration=1:timesteps # iterate for ts time step
      new_part = create_particles(source, iteration, 3)
      add!(new_part,part)
      
      ρ = particle_to_grid(part, grid)
      Ex, Ey = calculate_electric_field(solver, ρ, iteration)
      pE = grid_to_particle(grid, part, Ex, Ey)

      push_particles!(pusher, part, pE, Δt)

      # remove particles outside the computational domain
      p = 1
      while p <= part.np
        i, j, hx, hy = particle_cell(part.x[1:part.np,:], p, Δh)
        if i < 1 || i >= nx || j < 1 || j >= ny
          remove!(part, p)
          continue                           # repeat for the new particle
        end
        p = p + 1
      end

      after_loop(iteration)

      println(@sprintf "Time Step #%d, Particles #%d" iteration part.np)
    end

    exit_loop()
    println("Complete!")
  end
end
