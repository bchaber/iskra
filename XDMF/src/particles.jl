function save_species(species, it, fname,temporal)
  particles = new_child(temporal, "Grid")
  time      = new_child(particles, "Time")
  topology  = new_child(particles, "Topology")
  geometry  = new_child(particles, "Geometry")

  set_attributes(time, Value=h5attr(it, "time"))

  pt = it["particles"]
  for n in species
    global np = add_species(particles, fname, n, pt[n], geometry)
    set_attributes(particles, Name=n*" Particles", GridType="Uniform")
    set_attributes(time, Value=h5attr(it, "time"))
    set_attributes(topology, TopologyType="Polyvertex",
      NodesPerElement="1", NumberOfElements=np)
    set_attributes(geometry, GeometryType="X_Y_Z")
    break
  end
end

function add_species(particles, fname::String, n::String, g::HDF5.HDF5Group, geometry)
  np = length(g["id"])

  for m in HDF5.names(g["position"])
    d = g["position/"*m]
    dataitem = new_child(geometry, "DataItem")
    set_attributes(dataitem, Name=m, Format="HDF5", NumberType="Float", Precision="8", Dimensions=np)
    add_text(dataitem, @sprintf "%s:%s" fname HDF5.name(d))
  end

  for m in ("id",)
    d = g[m]
    attribute = new_child(particles, "Attribute")
    set_attributes(attribute, Name=n*m, AttributeType="Scalar", Center="Node")
    dataitem = new_child(attribute, "DataItem")
    set_attributes(dataitem, Name=m, Format="HDF5", NumberType="UInt", Dimensions=np)
    add_text(dataitem, @sprintf "%s:%s" fname HDF5.name(d))
  end
  return np
end