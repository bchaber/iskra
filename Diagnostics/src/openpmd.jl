ind(d::Unitful.Dimension{:Length})      = 1
ind(d::Unitful.Dimension{:Mass})        = 2
ind(d::Unitful.Dimension{:Time})        = 3
ind(d::Unitful.Dimension{:Current})     = 4
ind(d::Unitful.Dimension{:Temperature}) = 5
ind(d::Unitful.Dimension{:Amount})      = 6
ind(d::Unitful.Dimension{:Luminosity})  = 7

geo(grid::CartesianGrid{2}) = ("x", "y"), "cartesian", 2
geo(grid::CartesianGrid{3}) = ("x", "y", "z"), "cartesian", 3
geo(grid::AxialGrid{1})     = ("r",), "spherical", 2
geo(grid::AxialGrid{2})     = ("r", "z"), "cylindrical", 2
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
end

struct ParticleRecord{T, D} <: Record
  name       :: String
  input      :: Array{T, D}
  components :: Dict{Char, AbstractArray}
  metadata   :: ParticleRecordMetadata
end

function recreate(val::PaddedView, species)
  PaddedView(val.fillvalue, parent(val), (1:species.np,))
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

struct FieldRecord{T, D, N} <: Record
  input      :: Array{T, D}
  components :: Dict{Char, AbstractArray}
  metadata   :: FieldRecordMetadata{N}
end
function update!(record::FieldRecord, units, data; grid, optional...)
  record.input .= data
end
function zview(parent, indices)
  PaddedView(parent |> eltype |> zero, parent, (indices,))
end

function init!(::Type{ParticleRecord}, data, units;
  species, offset=0.0, weighted=false, withcomponents=false)
  unitDimension, unitSI = usi(units)

  name = species.name
  weightingPower = 1.0
  input = zero(data)

  zz = data |> eltype |> zero
  components = Dict{Char, AbstractArray}()
  if withcomponents
    components['x'] = view(input, 1:species.np, 1)
    components['y'] = view(input, 1:species.np, 2)
    components['z'] = view(input, 1:species.np, 3)
  elseif length(data) > 1
    components[' ']  = view(input, 1:species.np, 1)
  else
    components[' ']  = PaddedView(zz, input, (1:species.np,))
  end

  metadata = ParticleRecordMetadata(unitDimension, offset,
    weighted, weightingPower, unitSI)
  record =   ParticleRecord(name, input, components, metadata)
end

function init!(::Type{FieldRecord}, data, units;
  grid, offset=0.0, pos=nothing, withcomponents=false)
  unitDimension, unitSI   = usi(units)
  
  axisLabels, geometry, N = geo(grid)
  gridSpacing = grid.Δh[1:N]
  gridGlobalOffset = tuple(zeros(N)...)
  input = zero(data)

  zz = data |> eltype |> zero
  components = Dict{Char, AbstractArray}()
  if withcomponents
    components['x'] = view(input, :, :, 1)
    components['y'] = view(input, :, :, 2)
   #components['z'] = view(input, :, :, 3)
  elseif length(data) > 1
    components[' '] = view(input, :, :, 1)
  end

  if isnothing(pos)
    pos = tuple(zeros(N)...)
  end

  metadata = FieldRecordMetadata(unitDimension, offset,
    axisLabels, "C", geometry, "",
    gridGlobalOffset, gridSpacing, 1.0,
    "none",
    pos, unitSI)
  record   = FieldRecord(input, components, metadata)
end
