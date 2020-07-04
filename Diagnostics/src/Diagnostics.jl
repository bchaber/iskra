module Diagnostics

using Printf
using HDF5
using Unitful
using RegularGrid

import Dates

export @particle, @field
export @probe
export new_iteration, save_diagnostic

abstract type Record end

include("openpmd.jl")
include("circuit.jl")

const records = Dict{String, Record}()
const root = RootMetadata()
const fields = FieldMetadata()
const particles = ParticleMetadata()

function register_diagnostic(key::String, units::String, data, record; kwargs...)
  if haskey(records, key) == false
    println("Registering diagnostics ", key)
    records[key] = init!(record, data, units; kwargs...)
  end
  update!(records[key], units, data; kwargs...)
  nothing
end

macro field(key, units, data, grid, optional...)
  kwargs = [esc(arg) for arg in optional]
  quote
    Diagnostics.register_diagnostic($(esc(key)), $(esc(units)), $(esc(data)),
      FieldRecord, grid=$(esc(grid)); $(kwargs...))
  end
end

macro particle(key, units, data, part, optional...)
  kwargs = [esc(arg) for arg in optional]
  quote
    Diagnostics.register_diagnostic($(esc(key)), $(esc(units)), $(esc(data)),
      ParticleRecord, species=$(esc(part)); $(kwargs...))
  end
end

function forhdf5(v) v end
function forhdf5(v::NTuple{N,T}) where {N,T}
  N > 0 ? [v...] : Float64[]
end

function addattribute(sym, node, field)
  attr = attrs(node)
  attr[string(sym)] = forhdf5(field)
end

function addattributes(metadata, node; except=(), fields=nothing)
  attr = attrs(node)

  if isnothing(fields)
    fields = metadata |> typeof |> fieldnames
  end

  for sym in fields
    if sym ∈ except
      continue
    end
    field = getfield(metadata, sym)
    attr[string(sym)] = forhdf5(field)
  end
end

function new_iteration(f::Function, prefix, i, t, dt, directory="/tmp")
  iteration = Iteration(dt, t, 1.0)

  mkpath(@sprintf "%s/%s/hdf5" directory prefix)
  rootnode = h5open((@sprintf "%s/%s/hdf5/data%d.h5" directory prefix i), "w")
  basepath = g_create(rootnode, (@sprintf "data/%d" i))
  meshesnode = g_create(basepath, root.meshesPath)
  particlesnode = g_create(basepath, root.particlesPath)
  
  addattributes(root, rootnode)
  addattributes(fields, meshesnode)
  addattributes(iteration, basepath)
  addattributes(particles, particlesnode)

  try
    f(basepath)
  finally
    close(rootnode)
  end
end

function save_record(it, key, record)
  println("Cannot save abstract diagnostic data")
end

function save_record(it::HDF5Group, key::String, record::ParticleRecord)
  f = root.particlesPath * key

  for (component, output) in record.components
    g = component == ' ' ? f : f * '/' * component
    if length(record.input) > 1
      it[g] = output
    else
      g_create(it, f)
      addattribute(:value, it[g], record.input[1])
      addattribute(:shape, it[g], output |> size)
    end
    addattributes(record.metadata, it[g]; fields=(:unitSI,))
  end
  addattributes(record.metadata, it[f]; except=(:unitSI,))
end

function save_record(it::HDF5Group, key::String, record::FieldRecord)
  f = root.meshesPath * key

  for (component, output) in record.components
    g = component == ' ' ? f : f * '/' * component
    it[g] = output
    addattributes(record.metadata, it[g]; fields=(:position, :unitSI))
  end
  addattributes(record.metadata, it[f]; except=(:position, :unitSI))
end

function save_records(it, prefix)
  for key in keys(records)
    if startswith(key, prefix)
      save_record(it, key, records[key])
    end
  end
end

function save_record(it, key)
  if haskey(records, key)
    save_record(it, key, records[key])
  else
    println("Couldn't find diagnostic ", key)
  end
end

end
