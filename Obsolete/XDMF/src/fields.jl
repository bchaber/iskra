function write_fields(xdmf, xdoc)
  iteration = xdmf.iteration
  fname = xdmf.file.filename

  xroot = xdoc |> root
  domain = xroot["Domain"] |> first
  temporal = domain["Grid"] |> first
  
  fields = new_child(temporal, "Grid")
  time   = new_child(fields, "Time")
  topology  = new_child(fields, "Topology")
  geometry  = new_child(fields, "Geometry")

  set_attributes(fields, Name="Fields", GridType="Uniform")
  set_attributes(topology, TopologyType="3DCoRectMesh")
  set_attributes(geometry, GeometryType="ORIGIN_DXDYDZ")
  origin = new_child(geometry, "DataItem")
  spacing = new_child(geometry, "DataItem")
  set_attributes(origin,  Name="Origin",  Dimensions="3", NumberType="Float", Precision="4", Format="XML")
  set_attributes(spacing, Name="Spacing", Dimensions="3", NumberType="Float", Precision="4", Format="XML")

  it = xdmf.file[@sprintf "/data/%d" iteration]
  fl = xdmf.file[@sprintf "/data/%d/fields" iteration]
  for n in keys(fl)
    o = h5attr(fl[n], "gridGlobalOffset")
    s = h5attr(fl[n], "gridSpacing")
    global o, s, d = add_field(fields, fname, n, fl[n], o, s)
  end

  add_text(origin, o)
  add_text(spacing, s)
  set_attributes(topology, Dimensions=d)
  set_attributes(time, Value=h5attr(it, "time"))
end

function write_probes(xdmf, xdoc)
  iteration = xdmf.iteration
  fname = xdmf.file.filename

  xroot = xdoc |> root
  domain = xroot["Domain"] |> first
  temporal = domain["Grid"] |> first
  
  fields = new_child(temporal, "Grid")
  time   = new_child(fields, "Time")
  topology  = new_child(fields, "Topology")
  geometry  = new_child(fields, "Geometry")

  origin = new_child(geometry, "DataItem")
  spacing = new_child(geometry, "DataItem")

  set_attributes(fields, Name="Fields", GridType="Uniform")
  set_attributes(origin,  Name="Origin",  Dimensions="1", NumberType="Float", Precision="4", Format="XML")
  set_attributes(spacing, Name="Spacing", Dimensions="1", NumberType="Float", Precision="4", Format="XML")
  set_attributes(topology, TopologyType="Polyvertex")
  set_attributes(geometry, GeometryType="X_Y_Z")
  
  it = xdmf.file[@sprintf "/data/%d" iteration]
  fl = xdmf.file[@sprintf "/data/%d/fields" iteration]
  
  set_attributes(topology, Dimensions="1")
  set_attributes(time, Value=h5attr(it, "time"))

  for n in keys(fl)
      attribute = new_child(fields, "Attribute")
      set_attributes(attribute;
        Name=n, AttributeType="Scalar", Center="Node")
      dataitem = new_child(attribute, "DataItem")
      set_attributes(dataitem;
        Format="HDF5", NumberType="Float", Precision="8",
        Dimensions="1")
      add_text(dataitem, @sprintf "%s/%s:%s" pwd() fname name(fl[n]))
      #add_text(dataitem, "3.14")
  end
  for m in ("x","y","z")
    dataitem = new_child(geometry, "DataItem")
    set_attributes(dataitem, Name=m, Format="XML", NumberType="Float", Precision="8", Dimensions="1")
    add_text(dataitem, "0")
  end
  add_text(origin,  "0")
  add_text(spacing, "0")
end

function add_field(fields, fname::String, n::String, g::Dataset, origin, spacing)
  dx, dy = size(g)
  o = @sprintf "0.0 %g %g" origin...
  s = @sprintf "0.0 %g %g" spacing...
  d = @sprintf "1 %d %d" dy dx
  attribute = new_child(fields, "Attribute")
  set_attributes(attribute;
    Name=n, AttributeType="Scalar", Center="Node")
  dataitem = new_child(attribute, "DataItem")
  set_attributes(dataitem;
    Format="HDF5", NumberType="Float", Precision="8",
    Dimensions=d)
  add_text(dataitem, @sprintf "%s/%s:%s" pwd() fname name(g))
  return o, s, d
end

function add_field(fields, fname::String, n::String, g::Group,
  origin, spacing)
  for m in keys(g)
    global o, s, d = add_field(fields, fname, n*m, g[m], origin, spacing)
  end
  return o, s, d
end
