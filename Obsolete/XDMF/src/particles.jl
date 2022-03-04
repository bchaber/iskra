function write_species(xdmf, xdoc, species)
  iteration = xdmf.iteration
  fname = xdmf.file.filename

  xroot = xdoc |> root
  domain = xroot["Domain"] |> first
  temporal = domain["Grid"] |> first

  particles = new_child(temporal, "Grid")
  time      = new_child(particles, "Time")
  topology  = new_child(particles, "Topology")
  geometry  = new_child(particles, "Geometry")
  set_attributes(particles, Name=species*" Particles", GridType="Uniform")

  it = xdmf.file[@sprintf "/data/%d" iteration]
  pt = xdmf.file[@sprintf "/data/%d/particles" iteration]
  np = add_species(particles, fname, species, pt[species], geometry)
  set_attributes(time, Value=h5attr(it, "time"))
  set_attributes(topology, TopologyType="Polyvertex",
    NodesPerElement="1", NumberOfElements=np)
  set_attributes(geometry, GeometryType="X_Y_Z")
end

function add_species(particles, fname::String, n::String, g::Group, geometry)
  np = length(read(g["id"]))

  for m in keys(g["position"])
    d = g["position/"*m]
    dataitem = new_child(geometry, "DataItem")
    set_attributes(dataitem, Name=m, Format="HDF5", NumberType="Float", Precision="8", Dimensions=np)
    add_text(dataitem, @sprintf "%s/%s:%s" pwd() fname name(d))
  end

  for m in ("id",)
    d = g[m]
    attribute = new_child(particles, "Attribute")
    set_attributes(attribute, Name=n*m, AttributeType="Scalar", Center="Node")
    dataitem = new_child(attribute, "DataItem")
    set_attributes(dataitem, Name=m, Format="HDF5", NumberType="UInt", Dimensions=np)
    add_text(dataitem, @sprintf "%s/%s:%s" pwd() fname name(d))
  end
  for m in ("momentum/x","momentum/y","momentum/z")
    d = g[m]
    attribute = new_child(particles, "Attribute")
    set_attributes(attribute, Name=n*m, AttributeType="Scalar", Center="Node")
    dataitem = new_child(attribute, "DataItem")
    set_attributes(dataitem, Name=m, Format="HDF5", NumberType="Float", Precision="8", Dimensions=np)
    add_text(dataitem, @sprintf "%s/%s:%s" pwd() fname name(d))
  end
  return np
end
