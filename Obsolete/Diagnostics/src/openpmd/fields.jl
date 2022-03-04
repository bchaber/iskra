struct FieldMetadata
  fieldSolver :: String
  fieldSolverParameters :: String
  fieldBoundary :: NTuple{4, String}
  particleBoundary :: NTuple{4, String}
  currentSmoothing :: String
  chargeCorrection :: String
  FieldMetadata() = new("other", "Nagel",
    ("open","open","open","open"),
    ("periodic","periodic","periodic","periodic"),
    "none", "none")
end

struct FieldRecordMetadata{N}
  unitDimension :: NTuple{7, Float64}
  timeOffset :: Float64

  axisLabels :: NTuple{N, Char}
  dataOrder :: String
  geometry :: String
  geometryParameters :: String
  gridGlobalOffset :: NTuple{N, Float64}
  gridSpacing :: NTuple{N, Float64}
  gridUnitSI :: Float64

  fieldSmoothing :: String

  position :: NTuple{N, Float64}
  unitSI :: Float64
end

struct FieldRecord{T, D, N} <: Record
  input      :: Array{T, D}
  components :: Dict{Char, AbstractArray}
  metadata   :: FieldRecordMetadata{N}
end

function init!(::Type{FieldRecord}, data, units;
  grid, offset=0.0, pos=nothing, withcomponents=false)
  unitDimension, unitSI   = usi(units)
  
  axisLabels, geometry, L = geo(grid)
  gridSpacing = grid.Δh
  gridGlobalOffset = tuple(zeros(L)...)
  input = zero(data)

  components = Dict{Char, AbstractArray}()
  if withcomponents
    if ndims(data) == 2 # HACK!
      N, K = size(data)
      input = zeros(N, 1, K)
    end
    N, M, K = size(input)
    zz = zeros(N, M)
    components['x'] = K > 0 ? view(input, :, :, 1) : zz
    components['y'] = K > 1 ? view(input, :, :, 2) : zz
    components['z'] = K > 2 ? view(input, :, :, 3) : zz
  elseif length(data) > 1
    components[' '] = view(input, :, :, 1)
  end

  if isnothing(pos)
    pos = tuple(zeros(L)...)
  end

  metadata = FieldRecordMetadata(unitDimension, offset,
    axisLabels, "C", geometry, "",
    gridGlobalOffset, gridSpacing, 1.0,
    "none",
    pos, unitSI)
  record   = FieldRecord(input, components, metadata)
end

function update!(record::FieldRecord, units, data; grid, optional...)
  record.input[:] .= data[:] # HACK!
end