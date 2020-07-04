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