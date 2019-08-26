module ParticleInCell
  using FiniteDifferenceMethod
  
  include("sugar.jl")

  include("pic/macroparticle.jl")
  include("pic/cloud_in_cell.jl")
  include("pic/diagnostics.jl")
  include("pic/pushers.jl")
  include("pic/sources.jl")

  function remove_particles!(part, Δh, matches)
    p = 1
    px = view(part.x, 1:part.np, :)
    while p <= part.np
      i, j, _, _ = particle_cell(px, p, Δh)
      if matches(i,j)
        remove!(part, p)
        continue
      end
      p = p + 1
    end
  end

  function particle_cell(px, p, Δh)
    fij = @views 1 .+ px[p, :] ./ Δh
    ij = floor.(Int64, fij)
    hxy = fij .- ij
    i, j  = ij
    hx,hy = hxy
    return i, j, hx, hy
  end

  function cell_to_grid(cells, grid)
    nx, ny = size(grid)
    v = cells[1:nx,1:ny  ] .+ cells[2:nx+1,1:ny  ] .+
        cells[1:nx,2:ny+1] .+ cells[2:nx+1,2:ny+1]
    v ./= 4
  end

  # hooks
  function enter_loop() end
  function after_loop(it) end
  function exit_loop() end

  function init(src, Δt)
    part = src.species
    px, pv = @views part.x[1+part.np:end,:], part.v[1+part.np:end,:]
    n = create_particles!(src, px, pv, Δt)
    part.np += n
  end

  function solve(config, Δt=1e-5, timesteps=200)
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
    ε = cell_to_grid(cells["ε"], grid)
    for iteration=1:timesteps # iterate for ts time step
      # Create particles
      for src in sources
        part = src.species
        px, pv = @views part.x[1+part.np:end,:], part.v[1+part.np:end,:]
        np = create_particles!(src, px, pv, Δt)
        part.np += np
      end
      # Remove particles
      for part in species
        partE = grid_to_particle(grid, part, (i,j) -> E[i, j, :])
        push_particles!(pusher, part, partE, Δt)
        remove_particles!(part, Δh, (i,j) -> i < 1 || i >= nx || j < 1 || j >= ny)

        @diag "pv"*part.name ParticleVectorData(part.x,part.v,part.id, part.np)
        @diag "pE"*part.name ParticleVectorData(part.x,partE, part.id, part.np)
      end
      # Calculate charge density
      fill!(ρ, 0.0)
      for part in species
        n = particle_to_grid(part, grid, (p) -> part.np2c)
        ρ .+= n .* part.q

        @diag "n"*part.name NodeData(n, origin, spacing)
      end
      # Calculate electric field
      ϕ  = calculate_electric_potential(solver, ρ./ε)
      E  = calculate_electric_field(solver, ϕ)

      @diag "ρ" NodeData(ρ, origin, spacing)
      @diag "ϕ" NodeData(ϕ, origin, spacing)
      @diag "E" GridData(E, grid.x,  grid.y)

      after_loop(iteration)

      println("Time Step #", iteration, ", Particles #", [part.np for part in config.species])
    end

    exit_loop()
    println("Complete!")
  end
end
