mutable struct CircuitProbeRecord <: Record
  data  :: Matrix{Float64}
  metadata :: FieldRecordMetadata{2}
end

function CircuitProbeRecord(data, units::String; offset=0.0)
  unitDimension, unitSI   = usi(units)
  
  axisLabels, geometry, N = ("",""), "cartesian", 2
  gridSpacing = (.0,.0)
  gridGlobalOffset = (.0,.0)
  pos = (.0,.0)
  metadata = FieldRecordMetadata(unitDimension, offset,
    axisLabels, "C", geometry, "",
    gridGlobalOffset, gridSpacing, 1.0,
    "none",
    pos, unitSI)
  CircuitProbeRecord(data * ones(1,1), metadata)
end

function update!(record::CircuitProbeRecord, units, data; optional...)
  record.data .= data
end

macro probe(key, data, units, optional...)
  kwargs = [esc(arg) for arg in optional]
  quote
    Diagnostics.register_diagnostic($(esc(key)), $(esc(units)), $(esc(data)),
      CircuitProbeRecord; $(kwargs...))
  end
end

function save_record(it::HDF5Group, key::String, record::CircuitProbeRecord)
  f = @sprintf "%s/%s" root.meshesPath key
  it[f] = record.data
  addattributes(record.metadata, it[f])
end