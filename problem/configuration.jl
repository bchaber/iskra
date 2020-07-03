using ParticleInCell
using FiniteDifferenceMethod
using RegularGrid

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
  i,j,k = find_reference_node(nodes)

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

function create_gamma_ionization_source(rate, x, v)
  MaxwellianSource(rate, x, v)
end

function create_kinetic_species(name, N, q, m, weight)
  species = KineticSpecies{2,3}(name, N)
  species.q = q
  species.m = m
  species.wg .*= weight
  species.w0 = weight
  return species
end

function create_fluid_species(name, μ, q, m, mx, my)
  T = 300 # K
  n = zeros(mx, my)
  species = FluidSpecies(name, μ, q, m, n, T)
  return species
end

Config() = Config(nothing, nothing, nothing, [], [], [], nothing, nothing, nothing)
