import PlotVTK:scatter,heatmap,plot
import PlotVTK:pvd_create,pvd_save,pvd_add_timestep

macro diagnostics_init()
  quote
    global pvd_scatter = pvd_create("/tmp/scatter")
    global pvd_heatmap = pvd_create("/tmp/heatmap")
  end
end

macro diagnostics_save()
  quote
    pvd_save(pvd_scatter)
    pvd_save(pvd_heatmap)
  end
end

macro diagnostics_scatter(it, pE, part)
  quote
    local _part = $(esc(part))
    local _pE = $(esc(pE))
    local _it = $(esc(it)) 

    pxy = _part.x[1:_part.np,:]
    vxy = _part.v[1:_part.np,:]
    px, py, pz  = pxy[:,1], pxy[:,2], zeros(_part.np)
    vx, vy, vz  = vxy[:,1], vxy[:,2], zeros(_part.np)
    pEx,pEy,pEz =_pE[:,1], _pE[:,2],  zeros(_part.np)
    pvd_add_timestep(pvd_scatter,
      scatter(px, py, "/tmp/pE",
        "v" => (vx,vy,vz),
        "E" => (pEx,pEy,pEz);
      it=_it, save=false),
    _it)
  end
end

macro diagnostics_heatmap(it, spacing, origin, Ex, Ey, ϕ, ρ)
  quote
    local _Ex = $(esc(Ex)) 
    local _Ey = $(esc(Ey)) 
    local _ϕ  = $(esc(ϕ)) 
    local _ρ  = $(esc(ρ)) 
    local _it = $(esc(it)) 
    local _spacing = $(esc(spacing)) 
    local _origin  = $(esc(origin)) 

    pvd_add_timestep(pvd_heatmap, heatmap("Ex" =>_Ex, "/tmp/Ex", spacing=_spacing, origin=_origin; it=_it, save=false), _it)
    pvd_add_timestep(pvd_heatmap, heatmap("Ey" =>_Ey, "/tmp/Ey", spacing=_spacing, origin=_origin; it=_it, save=false), _it)
    pvd_add_timestep(pvd_heatmap, heatmap("ϕ" =>_ϕ, "/tmp/phi",  spacing=_spacing, origin=_origin; it=_it, save=false), _it)
    pvd_add_timestep(pvd_heatmap, heatmap("ρ" =>_ρ, "/tmp/rho",  spacing=_spacing, origin=_origin; it=_it, save=false), _it)
  end
end
