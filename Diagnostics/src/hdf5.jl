const root = RootMetadata()
const fields = FieldMetadata()
const particles = ParticleMetadata()

function HDF5.write_attribute(attr::HDF5.Attribute, memtype::HDF5.Datatype, strs::NTuple{N,Char}) where N
  HDF5.write_attribute(attr, memtype, join(strs))
end

function HDF5.dataspace(v::NTuple{N,Char}) where N
     HDF5.Dataspace(HDF5.h5s_create_simple(1, HDF5.hsize_t[N], HDF5.hsize_t[N]))
end

function HDF5.datatype(::NTuple{N,Char}) where N
    type_id = HDF5.h5t_copy(HDF5.H5T_C_S1)
    HDF5.h5t_set_size(type_id, 1)
    HDF5.h5t_set_cset(type_id, HDF5.H5T_CSET_UTF8)
    HDF5.Datatype(type_id)
end

function Base.setindex!(attr::HDF5.Attributes, val::NTuple{N, Int64}, name::String) where N
  setindex!(attr, [val...], name)
end

function Base.setindex!(attr::HDF5.Attributes, val::NTuple{N, Float64}, name::String) where N
  setindex!(attr, [val...], name)
end

function Base.setindex!(attr::HDF5.Attributes, val::NTuple{N, String}, name::String) where N
  setindex!(attr, [val...], name)
end

function addattribute(sym, node, field)
  attr = attributes(node)
  attr[string(sym)] = field
end

function addattributes(metadata, node; except=(), fields=nothing)
  attr = attributes(node)

  if isnothing(fields)
    fields = metadata |> typeof |> fieldnames
  end

  for sym in fields
    if sym ∈ except
      continue
    end
    field = getfield(metadata, sym)
    attr[string(sym)] = field
  end
end

function new_iteration(f::Function, prefix, i, t, dt)
  iteration = Iteration(dt, t, 1.0)

  mkpath(prefix * "/hdf5")
  rootnode = h5open((@sprintf "%s/hdf5/data%d.h5" prefix i), "w")
  basepath = create_group(rootnode, (@sprintf "data/%d" i))
  meshesnode = create_group(basepath, root.meshesPath)
  particlesnode = create_group(basepath, root.particlesPath)
  
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


function save_record(it::HDF5.Group, key::String, record::ParticleRecord)
  f = root.particlesPath * key

  for (component, output) in record.components
    g = component == ' ' ? f : f * '/' * component
    if length(record.input) > 1
      it[g] = output
    else
      create_group(it, f)
      addattribute(:value, it[g], record.input[1])
      addattribute(:shape, it[g], output |> size)
    end
    addattributes(record.metadata, it[g]; fields=(:unitSI,))
  end
  addattributes(record.metadata, it[f]; except=(:unitSI,))
end

function save_record(it::HDF5.Group, key::String, record::FieldRecord)
  f = root.meshesPath * key

  for (component, output) in record.components
    g = component == ' ' ? f : f * '/' * component
    it[g] = output
    addattributes(record.metadata, it[g]; fields=(:position, :unitSI))
  end
  addattributes(record.metadata, it[f]; except=(:position, :unitSI))
end

function save_record(it::HDF5.Group, key::String, record::CircuitProbeRecord)
  f = @sprintf "%s/%s" root.meshesPath key
  it[f] = record.data
  addattributes(record.metadata, it[f])
end
