module ParticleInCell
  export solve

  import Printf:@sprintf,println

  include("sugar.jl")
  include("diagnostics.jl")

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

    heatmap("dof" => reshape(1:length(grid),nx,ny), "/tmp/grid"; origin=origin,spacing=spacing)

@diagnostics_init

    for iteration=1:timesteps # iterate for ts time step
      new_part = create_particles(source, iteration, 3)
      add!(new_part,part)
      
      ρ = particle_to_grid(part, grid)
      Ex, Ey = calculate_electric_field(solver, ρ)
      pE = grid_to_particle(grid, part, Ex, Ey)

@diagnostics_scatter iteration pE part
@diagnostics_heatmap iteration spacing origin Ex Ey ρ

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

      println(@sprintf "Time Step #%d, Particles #%d" iteration part.np)
    end

@diagnostics_save

    println("Complete!")
  end
end
