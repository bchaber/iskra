module XDMF

using Printf
using LightXML

import HDF5
export xdmf
export save_fields, save_species

h5attr(g, a) = read(HDF5.attrs(g)[a])
include("fields.jl")
include("particles.jl")

function xdmf(func::Function, input, iterations, output)
  xdoc = XMLDocument()
  xdmf = create_root(xdoc, "Xdmf") 
  domain = new_child(xdmf, "Domain")
  temporal = new_child(domain, "Grid")

  set_attributes(xdmf, Version="3.0")
  set_attributes(temporal, GridType="Collection", CollectionType="Temporal")
  
  for i in iterations
    fname = @sprintf "%s%d.h5" input i
    iname = @sprintf "/data/%d" i
    
    f = HDF5.h5open(fname, "r");
    try
      func(f[iname], fname, temporal)
    finally
      close(f)
    end
  end

  save_file(xdoc, output)
end

end





