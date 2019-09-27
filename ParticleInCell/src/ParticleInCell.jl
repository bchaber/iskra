module ParticleInCell
  using LinearAlgebra
  using FiniteDifferenceMethod

  include("sugar.jl")

  include("pic/kinetic.jl")
  include("pic/fluid.jl")
  include("pic/cloud_in_cell.jl")
  include("pic/diagnostics.jl")
  include("pic/pushers.jl")
  include("pic/sources.jl")

  using Chemistry

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

  function init(src, species, Δt)
    sample!(src, species, Δt)
  end

  function advance!(part :: KineticSpecies, E, Δt, config)
    pusher = config.pusher
    grid = config.grid
    nx, ny = size(grid)

    partE = grid_to_particle(grid, part, (i,j) -> E[i, j, :])
    push_particles!(pusher, part, partE, Δt)
    remove_particles!(part, grid.Δh, (i,j) -> i < 1 || i >= nx || j < 1 || j >= ny)

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

  function maxwellian_velocity(ν)
    v = [randn(), randn(), 0]
    v ./ norm(v) .* ν
  end

  function perform!(collision::ElasticCollision, p, Δt, grid)
    ν = norm(collision.source.v[p,:])
    collision.source.v[p,:] .= maxwellian_velocity(ν)
  end

  function isotropic_velocity(ν)
    θ = 2rand() * π
    r = 2rand() - 1
    a = sqrt(1 - r^2)
    [cos(θ)*a, sin(θ)*a, r] .* ν
  end

  function thermal_velocity(T, m)
    kB = 1.3806503e-23
    sqrt(2kB*T/m)
  end

  function perform!(collision::IonizationCollision, p, Δt, grid)
    source = collision.source
    qe = 1.60217646e-19
    sv = view(source.v, p, :)
    sE = 0.5source.m*dot(sv, sv)/qe - 12.0697 #target.ionization_energy;
    if sE < 0
      println("Source species has not enough energy for ionization: ", sE)
      return
    end
    # randomly redistribute the remaining energy to the two electrons
    e1E = sE * rand()
    e2E = sE - e1E
    # speed reduced by the ionization energy
    e1ν = sqrt(e1E*qe*2/source.m)
    e2ν = sqrt(e2E*qe*2/source.m)
            
    sv .= isotropic_velocity(e1ν)

    # assume the new electron and ion are created at the neutral temperature
    T = 300 # target.temperature # 300K = 25°C
    # create new ion and electron
    source.x[source.np+1,:] .= source.x[p,:]
    source.v[source.np+1,:] .= isotropic_velocity(e2ν);
    source.np += 1
    for product in collision.products
      if product == source
        continue
      end
      mp = source.w0/product.w0 + rand()
      νth = thermal_velocity(T, product.m)
      for i=1:round(Integer, mp)
        product.x[product.np+1,:] .= source.x[p,:]
        product.v[product.np+1,:] .= maxwellian_velocity(νth);
        product.np += 1
      end
    end
  end

  function perform!(dsmc::DirectSimulationMonteCarlo, Δt, grid, E)
    ν = zeros(size(grid)) # collision count
    Δh = grid.Δh
    for collision in dsmc.collisions
      source, target = collision.source, collision.target
      cache(source, target, dsmc) # assign particles to cells
    end
    @diag "ν-mcc" NodeData(ν, grid.origin, [Δh,Δh])
  end

  function perform!(mcc::MonteCarloCollisions, Δt, grid, E)
    ν = zeros(size(grid)) # collision count
    Δh = grid.Δh
    for collision in mcc.collisions
      source, target = collision.source, collision.target
      for p=1:source.np
        i, j, _, _ = particle_cell(source.x, p, grid.Δh)
        n = density(target, grid)[i,j]
        if n < 0
          println("Density is negative, skipping")
          continue
        end
        sv = source.v
        tv = (target.q/target.m)*E*Δt
        gv = norm(tv[i,j,:] .- sv[p,:])
        σ  = collision.rate(gv)
        P  = @. 1 - exp(-σ * gv * Δt * n);
        R  = rand()
        if P < R
          continue
        end
        perform!(collision, p, Δt, grid)
        ν[i,j] += 1
      end
    end
    @diag "ν-mcc" NodeData(ν, grid.origin, [Δh,Δh])
  end

  function perform!(network::ChemicalReactionNetwork, Δt, grid, E)
    Δn = Dict()
    for reaction in network.reactions
      rate = reaction.rate.(E[:,:,1])
      for (species, coeff) in reaction.stoichiometry
        n = species.n
        if haskey(Δn, species) ≠ true
          Δn[species] = zeros(size(n))
        end
        concentration = ones(size(n))
        for (species, order) in reaction.reactants
          concentration .*= species.n.^order
        end
        Δn[species] .+= Δt .* coeff .* rate .* concentration
      end
    end

    for species in keys(Δn)
      species.n .+= Δn[species]
      @diag "Δn"*species.name NodeData(Δn[species], grid.origin, [1,1]*grid.Δh)
    end
  end

  function solve(config, Δt=1e-5, timesteps=200, ε0=1.0)
    pusher = config.pusher
    interactions = config.interactions
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
        part.n = density(part, grid)
        ρ .+= part.n .* part.q

        @diag "n"*part.name NodeData(part.n, origin, spacing)
      end
      # Solve reactions
      for interaction in interactions
        perform!(interaction, Δt, grid, E)
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
