import ParticleInCell
import FiniteDifferenceMethod

mutable struct Config
  solver
  pusher
  chemistry
  species
  sources
  grid
  cells
end

function set_permittivity(εr)
  grid = config.grid
  FiniteDifferenceMethod.create_generalized_poisson_solver(grid, εr)
end

function add_electrode(nodes, voltage)
  solver, dof = config.solver, config.grid.dof
  FiniteDifferenceMethod.apply_dirichlet(solver, dof[nodes], voltage)
end

function create_gamma_ionization_source(electrons, rate, x, v)
  ParticleInCell.create_maxwellian_source(electrons, rate, x, v)
end

function create_kinetic_species(name, N, q, m, np2c)
  species = ParticleInCell.KineticSpecies(name, N)
  species.q = q
  species.m = m
  species.np2c = np2c
  return species
end

function create_fluid_species(name, μ, q, m)
  species = ParticleInCell.FluidSpecies(name, μ, q, m, zeros(0, 0))
  return species
end

Config() = Config(nothing, nothing, nothing, [], [], nothing, nothing)
