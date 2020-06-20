function save_fields(it, fname, temporal)
  fields = new_child(temporal, "Grid")
  time   = new_child(fields, "Time")
  topology  = new_child(fields, "Topology")
  geometry  = new_child(fields, "Geometry")

  set_attributes(fields, Name="Fields", GridType="Uniform")
  set_attributes(time, Value=h5attr(it, "time"))
  set_attributes(topology, TopologyType="3DCoRectMesh")
  set_attributes(geometry, GeometryType="ORIGIN_DXDYDZ")

  origin = new_child(geometry, "DataItem")
  set_attributes(origin, Name="Origin", Dimensions="3", NumberType="Float", Precision="4", Format="XML")
  spacing = new_child(geometry, "DataItem")
  set_attributes(spacing, Name="Spacing", Dimensions="3", NumberType="Float", Precision="4", Format="XML")

  fl = it["fields"]
  for n in HDF5.names(fl)
    o  = h5attr(fl[n], "gridGlobalOffset")
    s = h5attr(fl[n], "gridSpacing")
    global o, s, d = add_field(fields, fname, n, fl[n], o, s)
  end

  add_text(origin, o)
  add_text(spacing, s)
  set_attributes(topology, Dimensions=d)
end

function add_field(fields, fname::String, n::String, g::HDF5.HDF5Dataset, origin, spacing)
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
  add_text(dataitem, @sprintf "%s:%s" fname HDF5.name(g))
  return o, s, d
end

function add_field(fields, fname::String, n::String, g::HDF5.HDF5Group,
  origin, spacing)
  for m in HDF5.names(g)
    global o, s, d = add_field(fields, fname, n*m, g[m], origin, spacing)
  end
  return o, s, d
end