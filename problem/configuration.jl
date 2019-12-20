import ParticleInCell
import FiniteDifferenceMethod

mutable struct Config
  solver
  pusher
  interactions
  species
  sources
  circuit
  tracker
  grid
  cells
end

function set_permittivity(εr)
  grid = config.grid
  FiniteDifferenceMethod.create_generalized_poisson_solver(grid, εr)
end

function add_electrode(nodes, voltage)
  solver = config.solver
  FiniteDifferenceMethod.apply_dirichlet(solver, nodes, voltage)
end

function create_gamma_ionization_source(rate, x, v)
  ParticleInCell.MaxwellianSource(rate, x, v)
end

function create_kinetic_species(name, N, q, m, weight)
  species = ParticleInCell.KineticSpecies(name, N)
  species.q = q
  species.m = m
  species.wg .*= weight
  species.w0 = weight
  return species
end

function create_fluid_species(name, μ, q, m, mx, my)
  T = 300 # K
  n = zeros(mx, my)
  species = ParticleInCell.FluidSpecies(name, μ, q, m, n, T)
  return species
end

Config() = Config(nothing, nothing, [], [], [], nothing, nothing, nothing, nothing)
