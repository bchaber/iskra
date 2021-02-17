using ParticleInCell
using FiniteDifferenceMethod
using RegularGrids

mutable struct Config
  solver
  pusher
  tracker
  interactions
  species
  sources
  circuit
  grid
  cells
end

function set_permittivity(εr)
  grid = config.grid
  FiniteDifferenceMethod.create_generalized_poisson_solver(grid, εr, ε0)
end

function create_electrode(nodes, config;   fixed=false, σ=0.0, ϕ=0.0)
  if config.grid == nothing
    println("No grid defined. EXIT")
  end

  if config.solver == nothing
    println("No field solver defined. EXIT")
  end

  if config.tracker == nothing
    config.tracker = ParticleInCell.create_surface_tracker(config.grid)
  end

  electrode = create_electrode(nodes, config.solver, config.grid;
    fixed=fixed, σ=σ, ϕ=ϕ)
  ParticleInCell.track_surface!(config.tracker, nodes, electrode)
  return electrode
end

function create_electrode(nodes, ps, grid::CartesianGrid{2}; fixed=false, σ=0.0, ϕ=0.0)
  nx, ny = size(grid)
  Δx, Δy = grid.Δh
  Δz = 1.0
  function calculate_area(nodes)
    area = 0.
    for c in findall(nodes)
      i, j  = c.I
      area += (i < nx && nodes[i+1,j]) ? Δx*Δz : 0
      area += (j < ny && nodes[i,j+1]) ? Δy*Δz : 0
    end
    return area
  end
  function find_reference_node(nodes)
    c = findfirst(isequal(true), nodes)
    c ≠ nothing ? c.I : nothing
  end
  area = calculate_area(nodes)
  i,j  = find_reference_node(nodes)

  if fixed
    apply_dirichlet(ps, nodes, ϕ)
    ϕ0 = get_rhs(ps, :ϕ, i, j) .= ϕ
    return ParticleInCell.FixedPotentialElectrode(ϕ0, 0.0, area)
  else
    dof = FiniteDifferenceMethod.add_new_dof(ps, :σ)
    apply_neumann(ps, nodes, dof)
    σ0 = get_rhs(ps, :σ, dof) .= σ
    ϕ0 = get_solution(ps, :ϕ, i, j)
    return ParticleInCell.FloatingPotentialElectrode(ϕ0, σ0, 0.0, area)
  end
end

function plasma_frequency(species::KineticSpecies, n)
  q, m = species.q, species.q
  2π*q*sqrt(n/m)
end

function thermal_speed(T, m)
  sqrt(2kB*T/m)
end

function create_gamma_ionization_source(species::KineticSpecies{D,V}, x; dx=nothing, T=300K, rate=1) where {D,V}
  vth = thermal_speed(T, species.m) * ones(1, V)
  if isnothing(dx) dx = zero(x) end
  MaxwellianSource{D,V}(float.(rate), float.(x), float.(vth); dx=float.(dx))
end

function create_thermalized_beam(species::KineticSpecies{D,V}, x, vb; dx=nothing, T=300K, rate=1) where {D,V}
  vth = thermal_speed(T, species.m) * ones(1, V)
  if isnothing(dx) dx = zero(x) end
  MaxwellianSource{D,V}(float.(rate), float.(x), float.(vth); dx=float.(dx), dv=float.(vb))
end

function create_kinetic_species(name, N, q, m, weight; D=2, V=3)
  species = KineticSpecies{D,V}(name, N)
  species.q = q
  species.m = m
  species.wg .*= weight
  species.w0 = weight
  return species
end

function create_fluid_species(name, μ, q, m, mx, my; T=300.0K)
  n = zeros(mx, my)
  species = FluidSpecies(name, μ, q, m, n, T)
  return species
end

function calculate_potential_energy(species)
  ts = 0
  Uk = zeros(ts)
  for part in species
    for i=1:ts
      mass = part.m
      vx = view(part.v, 1:part.np, 1)
      vy = view(part.v, 1:part.np, 2)
      vz = view(part.v, 1:part.np, 3)
      uk = vx.^2 .+ vy.^2 .+ vz.^2
      Uk[i] .+= 0.5sum(uk)*mass
    end    
  end
  return Uk
end

function calculate_kinetic_energy(E, grid)
  ts = 0
  Up = zeros(ts)
  volume = cell_volume(grid)
  for i=1:ts
    ε  = ε0
    Ex = view(E, :, :, 1)
    Ey = view(E, :, :, 2)
    Ez = view(E, :, :, 3)
    up = Ex.^2 .+ Ey.^2 .+ Ez.^2
    Up[i] .+= 0.5sum(volume .* up)*ε
  end
  return Up
end

Config() = Config(nothing, nothing, nothing, [], [], [], nothing, nothing, nothing)
