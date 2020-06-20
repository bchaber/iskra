ind(d::Unitful.Dimension{:Length})      = 1
ind(d::Unitful.Dimension{:Mass})        = 2
ind(d::Unitful.Dimension{:Time})        = 3
ind(d::Unitful.Dimension{:Current})     = 4
ind(d::Unitful.Dimension{:Temperature}) = 5
ind(d::Unitful.Dimension{:Amount})      = 6
ind(d::Unitful.Dimension{:Luminosity})  = 7

geo(grid::UniformGrid{XY2D}) = ("x", "y"), "cartesian", 2
geo(grid::UniformGrid{RZ2D}) = ("r", "z"), "cartesian", 2
current_version() = "git+?"
current_date() = Dates.format(Dates.now(), "Y/m/d HH:MM")

function spam(d::Unitful.Dimensions{N}) where {N}
  dimensions = zeros(Float64, 7)
  for d in N
    dimensions[ind(d)] = convert(Float64, d.power)
  end
  tuple(dimensions...)
end

function eggs(d)
  d
end

function eggs(d::Unitful.Quantity{T,D,U}) where {T,D,U}
  convert(Float64, d.val)
end

function usi(units::String)
  U = 1(units |> uparse)
  D = U |> dimension
  S = U |> upreferred
  spam(D), eggs(S)
end

struct Iteration
  dt :: Float64
  time :: Float64
  timeUnitSI :: Float64
end

struct RootMetadata
  openPMD :: String
  openPMDextension :: UInt32
  iterationEncoding :: String
  iterationFormat :: String
  basePath :: String
  meshesPath :: String
  particlesPath :: String
  author :: String
  software :: String
  softwareVersion :: String
  date :: String
end

RootMetadata() = RootMetadata("1.1.0", 1, "fileBased", "data%T.h5", "/data/%T",
  "fields/", "particles/", "Bartosz Chaber <bartosz.chaber@ee.pw.edu.pl>",
  "iskra", current_version(), current_date())

struct FieldMetadata
  fieldSolver :: String
  fieldSolverParameters :: String
  fieldBoundary :: NTuple{4, String}
  particleBoundary :: NTuple{4, String}
  currentSmoothing :: String
  chargeCorrection :: String
end
FieldMetadata() = FieldMetadata("other", "Nagel",
  ("open","open","open","open"),
  ("periodic","periodic","periodic","periodic"),
  "none", "none")

struct ParticleMetadata
  particleShape :: Float64
  currentDeposition :: String
  particlePush :: String
  particleInterpolation :: String
  particleSmoothing :: String
end
ParticleMetadata() = ParticleMetadata(0., "none", "Boris", "other", "none")

struct ParticleRecordMetadata
  unitDimension :: NTuple{7, Float64}
  timeOffset :: Float64

  macroWeighted :: UInt8
  weightingPower :: Float64

  unitSI :: Float64 # Component

  value :: Float64
  shape :: Int64
end

mutable struct ParticleRecord{T, N, M} <: Record
  data :: Array{T, N}
  components :: NTuple{M, String}
  name :: String
  np :: Int64
  metadata :: ParticleRecordMetadata
end
function update!(record::ParticleRecord, units, data; species, optional...)
  record.np    = species.np
  record.data .= data
end

struct FieldRecordMetadata{N}
  unitDimension :: NTuple{7, Float64}
  timeOffset :: Float64

  axisLabels :: NTuple{N, String}
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

struct FieldRecord{T, N, M, D} <: Record
  data  :: Array{T, D}
  components :: NTuple{M, String}
  metadata :: FieldRecordMetadata{N}
end
function update!(record::FieldRecord, units, data; grid, optional...)
  record.data .= data
end

@inline function checkfieldcomponents(data, components)
  N, M = size(data, 3), length(components)
  if N > 1 || M > 0
    @assert M == N "Number of components $M does not match data size $N"
  end
end

@inline function checkparticlecomponents(data, components)
  N, M = size(data, 2), length(components)
  if N > 1 || M > 0
    @assert M == N "Number of components $M does not match data size $N"
  end
end

function ParticleRecord(data, units;
  species, offset=0.0, weighted=false, components=())
  checkparticlecomponents(data, components)
  unitDimension, unitSI = usi(units)

  npar = species.np
  name = species.name
  weightingPower = 1.0

  value = length(data) > 1 ? NaN : data[1]
  shape = species.np
  
  metadata = ParticleRecordMetadata(unitDimension, offset,
    weighted, weightingPower, unitSI, value, shape)
  ParticleRecord(data, components, name, npar, metadata)
end

function FieldRecord(data, units;
  grid, offset=0.0, pos=nothing, components=())
  checkfieldcomponents(data, components)
  unitDimension, unitSI   = usi(units)
  
  axisLabels, geometry, N = geo(grid)
  gridSpacing = grid.Δh[1:N]
  gridGlobalOffset = tuple(zeros(N)...)

  if isnothing(pos)
    pos = tuple(zeros(N)...)
  end

  metadata = FieldRecordMetadata(unitDimension, offset,
    axisLabels, "C", geometry, "",
    gridGlobalOffset, gridSpacing, 1.0,
    "none",
    pos, unitSI)
  FieldRecord(data, components, metadata)
end
