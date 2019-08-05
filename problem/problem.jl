import Random
Random.seed!(0)
############################################
include("configuration.jl")
config = Config()
############################################
include("units_and_constants.jl")
const Mi = 99u
############################################
n0 = 1e12       # density in #/m^3
ϕ0 =-1V         # reference potential
ϕp =-1V         # wall potential
Te = 1.0eV      # electron temperature in eV
Ti = 0.1eV      # ion temperature in eV
v_drift = .5kmps  # ion injection velocity
# calculate plasma parameters
λD = sqrt(ɛ0*Te/(n0*e)) # Debye length
v_th = sqrt(2e*Ti/Mi)   # thermal velocity with Ti in eV
# set simulation domain
nx = 20             # number of nodes in x direction
ny = 20             # number of nodes in y direction
ts = 100            # number of time steps
Δh = 0.01m          # cell size
Δt = 0.001Δh/v_drift # time step, at vdrift move 0.10dx
nn = nx*ny          # total number of nodes
Lx = nx*Δh      # domain length in x direction
Ly = ny*Δh      # domain length in y direction
############################################
create_domain(0m:Δh:Lx, 0m:Δh:Ly)
bcs = zeros(Int8, size(config.grid))
inbox = (0.02m .<= config.grid.x .<= 0.05m) .&
        (0.02m .<= config.grid.y .<= 0.04m)
bcs[inbox].= 5;
bcs[ 1,10:ny] .= 1
bcs[nx, 1:10] .= 2
add_electrode(bcs .== 1, 1V)
add_electrode(bcs .== 2,-1V)
############################################
add_species("e-", 20_000, 1e, Mi, 100)
add_particle_source(3, Lx, Ly, v_th, v_drift)
############################################
import ParticleInCell: solve
@time solve(config)
