module Diagnostics

using Printf
using HDF5
using Unitful
using RegularGrids

import Dates

export @particle, @field
export @probe
export new_iteration, save_record, save_records

abstract type Record end

include("openpmd/root.jl")
include("openpmd/fields.jl")
include("openpmd/particles.jl")
include("circuit.jl")
include("hdf5.jl")

const records = Dict{String, Record}()

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

function save_record(it, key, record)
  println("Cannot save abstract diagnostic data")
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