import FiniteDifferenceMethod
import ParticleInCell
import RegularGrid

const RG = RegularGrid
const FDM = FiniteDifferenceMethod
const PIC = ParticleInCell
const LOAD = ParticleInCell.Loader

mutable struct Config
  field_solver
  loader
  ############
  cells
  grid
  species
end

function add_electrode(nodes, voltage)
  dof = config.grid.dof
  config.field_solver.apply_dirichlet(dof[nodes], voltage)
end

function add_species(name, N, q, m, np2c)
  species = PIC.Species(name, N)
  species.q = q
  species.m = m
  species.np2c = np2c
  config.species = species
end

function add_particle_source(max_np, xe, ye, vth, vdr)
  config.loader.setup(max_np, xe, ye, vth, vdr)
end

function create_domain(xs, ys)
  nx, ny = length(xs), length(ys)
  hx, hy = xs[2] - xs[1], ys[2] - ys[1]
  xm, ym = xs[2:nx] .- hx/2, ys[2:ny] .- hy/2

  config.cells = RG.UniformGrid(xs, ys)
  config.grid  = RG.UniformGrid(xm, ym)

  config.field_solver.setup(config.grid)
end

Config() = Config(FDM, LOAD, nothing, nothing, nothing)
