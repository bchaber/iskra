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

import Diagnostics
struct ParticleVectorData <: Diagnostics.DiagnosticData
  x :: Array{Float64,1}
  y :: Array{Float64,1}
  u :: Array{Float64,1}
  v :: Array{Float64,1}
  w :: Array{Float64,1}
 id :: Array{UInt64,1}
 it :: Integer
end

struct NodeData <: Diagnostics.DiagnosticData
  u :: Array{Float64,2}
 or :: Array{Float64,1}
 sp :: Array{Float64,1}
 it :: Integer
end

struct GridData <: Diagnostics.DiagnosticData
  u :: Array{Float64,3}
  x :: Array{Float64,2}
  y :: Array{Float64,2}
 it :: Integer
end

import PlotVTK: pvd_add_timestep, field_as_points, field_as_vectors, field_as_grid, field_as_vectors
Diagnostics.save_diagnostic(dname::String, d::ParticleVectorData, cname::String, c::Any, it::Integer) =
  pvd_add_timestep(c, field_as_vectors(d.x, d.y, dname, dname => (d.u, d.v, d.w), "uuid" => (d.id,), it=it, save=false), it)
Diagnostics.save_diagnostic(dname::String, d::NodeData, cname::String, c::Any, it::Integer) =
  pvd_add_timestep(c, field_as_points(dname  => d.u, dname, spacing=d.sp, origin=d.or, it=it, save=false), it)
Diagnostics.save_diagnostic(dname::String, d::GridData, cname::String, c::Any, it::Integer) =
  pvd_add_timestep(c, field_as_grid(d.x, d.y, dname  => d.u, dname, it=it, save=false), it)

  function remove_particles!(part, Δh, matches)
    p = 1
    while p <= part.np
      i, j, hx, hy = particle_cell(part.x[1:part.np,:], p, Δh)
      if matches(i,j)
        remove!(part, p)
        continue
      end
      p = p + 1
    end
  end

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
    sources = config.sources
    species = config.species
    solver = config.solver
    grid = config.grid

    nx, ny = size(grid)
    Δh = grid.Δh
    Δh² = Δh^2
    spacing = [1,1]*grid.Δh
    origin  = grid.origin

    enter_loop()

    for iteration=1:timesteps # iterate for ts time step
      ρ  = zeros(nx, ny)

      for src in sources
        create_particles!(src, iteration)
      end

      for part in species
        n   = particle_to_grid(part, grid, (p) -> part.np2c)
        ρ .+= particle_to_grid(part, grid, (p) -> part.np2c * part.q)
        Diagnostics.register_diagnostic("n"*part.name, NodeData(n, origin, spacing, iteration))
      end

      ϕ  = calculate_electric_potential(solver, ρ)
      E  = calculate_electric_field(solver, ϕ)

      Diagnostics.register_diagnostic("ρ", NodeData(ρ, origin, spacing, iteration))
      Diagnostics.register_diagnostic("ϕ", NodeData(ϕ, origin, spacing, iteration))
      Diagnostics.register_diagnostic("E", GridData(E, grid.x,  grid.y, iteration))

      for part in species
        partE = grid_to_particle(grid, part, (i,j) -> E[i, j, :])
        px, pv = part.x[1:part.np,:], part.v[1:part.np,:]
        pE = partE[1:part.np,:]
        pid = part.id[1:part.np]
        Diagnostics.register_diagnostic("pv"*part.name, ParticleVectorData(px[:,1],px[:,2],pv[:,1],pv[:,2],zeros(part.np),pid, 1))
        Diagnostics.register_diagnostic("pE"*part.name, ParticleVectorData(px[:,1],px[:,2],pE[:,1],pE[:,2],zeros(part.np),pid, 1))
        push_particles!(pusher, part, partE, Δt)
        remove_particles!(part, Δh, (i,j) -> i < 1 || i >= nx || j < 1 || j >= ny)
      end

      after_loop(iteration)

      println(@sprintf "Time Step #%d, Particles #%s" iteration [part.np for part in config.species])
    end

    exit_loop()
    println("Complete!")
  end
end
