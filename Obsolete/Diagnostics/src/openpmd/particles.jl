struct ParticleMetadata
  particleShape :: Float64
  currentDeposition :: String
  particlePush :: String
  particleInterpolation :: String
  particleSmoothing :: String
  ParticleMetadata() = new(0., "none", "Boris", "other", "none")
end

struct ParticleRecordMetadata
  unitDimension :: NTuple{7, Float64}
  timeOffset :: Float64

  macroWeighted :: UInt8
  weightingPower :: Float64

  unitSI :: Float64 # Component
end

struct ParticleRecord{T, D} <: Record
  name       :: String
  input      :: Array{T, D}
  components :: Dict{Char, AbstractArray}
  metadata   :: ParticleRecordMetadata
end

function init!(::Type{ParticleRecord}, data, units;
  species, offset=0.0, weighted=false, withcomponents=false)
  unitDimension, unitSI = usi(units)

  name = species.name
  weightingPower = 1.0
  input = zero(data)

  components = Dict{Char, AbstractArray}()
  if withcomponents
    N, M = size(input)
    zz = view(zeros(N, 1), 1:species.np, 1)
    components['x'] = M > 0 ? view(input, 1:species.np, 1) : zz
    components['y'] = M > 1 ? view(input, 1:species.np, 2) : zz
    components['z'] = M > 2 ? view(input, 1:species.np, 3) : zz
  elseif length(data) > 1
    components[' ']  = view(input, 1:species.np, 1)
  else
    components[' ']  = similar(input, species.np)
  end

  metadata = ParticleRecordMetadata(unitDimension, offset,
    weighted, weightingPower, unitSI)
  record =   ParticleRecord(name, input, components, metadata)
end

function recreate(val::Array, species)
  similar(val, species.np)
end

function recreate(val::SubArray, species)
  _, i = val.indices
  view(parent(val), 1:species.np, i)
end

function update!(record::ParticleRecord, units, data; species, optional...)
  record.input .= data
  for (key, val) in record.components
    record.components[key] = recreate(val, species)
  end
end