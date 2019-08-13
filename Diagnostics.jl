module Diagnostics

import Printf
import PlotVTK:pvd_create,pvd_save

abstract type DiagnosticData end

diagnostics = Dict{String, DiagnosticData}()
containers  = Dict{String, Any}()
directory = "/tmp/"

function open_container(cname::String)
  containers[cname] = pvd_create(directory*cname)
end

function register_diagnostic(dname::String, data::DiagnosticData)
  diagnostics[dname] = data
end

function save_diagnostic(dname::String, d::DiagnosticData, cname::String, container, iteration::Integer)
  Printf.println("Cannot save abstract diagnostic data")
end

function save_diagnostic(dname::String, cname::String, it::Integer)
  if haskey(diagnostics, dname)
    data = diagnostics[dname]
    container = containers[cname]
    save_diagnostic(directory*dname, data, directory*cname, container, it)
  else
    Printf.println("Couldn't find diagnostic "*dname)
  end
end

function close_container(cname::String)
  pvd_save(containers[cname])
  delete!(containers, cname)
end

end
