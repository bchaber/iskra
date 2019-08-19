import ParticleInCell
import FiniteDifferenceMethod

mutable struct Config
  solver
  pusher
  species
  sources
  grid
end

function add_electrode(nodes, voltage)
  solver, dof = config.solver, config.grid.dof
  FiniteDifferenceMethod.apply_dirichlet(solver, dof[nodes], voltage)
end

function create_gamma_ionization_source(electrons, rate, x, v)
  ParticleInCell.create_maxwellian_source(electrons, rate, x, v)
end

function create_species(name, N, q, m, np2c)
  species = ParticleInCell.Species(name, N)
  species.q = q
  species.m = m
  species.np2c = np2c
  return species
end

Config() = Config(nothing, nothing, nothing, nothing, nothing)
