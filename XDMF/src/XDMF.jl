module XDMF

using Printf

import LightXML: XMLDocument, XMLElement
import LightXML: root, create_root, new_child, set_attributes, add_text, find_element, save_file
import HDF5: HDF5File, HDF5Dataset, HDF5Group
import HDF5: attrs, h5open, name, names

export xdmf, new_document, save_document
export write_fields, write_species, write_probes

h5attr(g, a) = read(attrs(g)[a])

struct XDMFFile
  iteration :: Int64
  file :: HDF5File
end

include("fields.jl")
include("particles.jl")

function new_document()
  xdoc = XMLDocument()
  xroot = create_root(xdoc, "Xdmf") 
  domain = new_child(xroot, "Domain")
  temporal = new_child(domain, "Grid")
  set_attributes(xroot, Version="3.0")
  set_attributes(temporal, GridType="Collection", CollectionType="Temporal")
  return xdoc
end

function save_document(xdoc::XMLDocument, filename::String)
  save_file(xdoc, "./xdmf/"*filename*".xdmf")
end

function xdmf(func::Function, iterations)
  for i in iterations
    mkpath("./xdmf")
    name = @sprintf "./hdf5/data%d.h5" i
    file = h5open(name, "r");
    try
      XDMFFile(i, file) |> func
    finally
      close(file)
    end
  end
end

end