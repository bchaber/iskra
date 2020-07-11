const root = RootMetadata()
const fields = FieldMetadata()
const particles = ParticleMetadata()

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

function new_iteration(f::Function, prefix, i, t, dt)
  iteration = Iteration(dt, t, 1.0)

  mkpath(prefix * "/hdf5")
  rootnode = h5open((@sprintf "%s/hdf5/data%d.h5" prefix i), "w")
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

function save_record(it::HDF5Group, key::String, record::CircuitProbeRecord)
  f = @sprintf "%s/%s" root.meshesPath key
  it[f] = record.data
  addattributes(record.metadata, it[f])
end