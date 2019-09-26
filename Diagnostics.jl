module Diagnostics

import PlotVTK:pvd_create,pvd_save

abstract type DiagnosticData end

diagnostics = Dict{String, DiagnosticData}()
containers  = Dict{String, Any}()
directory = "/tmp/"

function open_container(cname::String)
  containers[cname] = pvd_create(directory*cname)
end

function register_diagnostic(dname::String, data::DiagnosticData)
  if !haskey(diagnostics, dname)
    println("New diagnostics registered: ", dname)
  end
  diagnostics[dname] = data
end

function save_diagnostic(dname::String, d::DiagnosticData, cname::String, container, iteration::Integer, t::Float64)
  println("Cannot save abstract diagnostic data")
end

function save_diagnostic(dname::String, cname::String, it::Integer, t=nothing)
  if haskey(diagnostics, dname)
    data = diagnostics[dname]
    container = containers[cname]
    if isnothing(t) t = it end
    save_diagnostic(dname, data, directory*cname*"_", container, it, float(t))
  else
    println("Couldn't find diagnostic "*dname)
  end
end

function close_container(cname::String)
  pvd_save(containers[cname])
  delete!(containers, cname)
end

macro diag(dname, data)
  quote
   Diagnostics.register_diagnostic($(esc(dname)), $(esc(data)))
  end
end

end
