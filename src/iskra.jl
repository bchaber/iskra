println("[*] Checking runtime dependencies")
import Pkg
Pkg.add("ProgressMeter")
Pkg.activate(".")
Pkg.resolve()

module IskraRuntime
  include("../problem/configuration.jl")
  include("../problem/units_and_constants.jl")

  using Random
  using ProgressMeter
  using ParticleInCell

  # define initial values
  grid = solver = pusher = nothing
  species = []
  interactions = []
  ParticleInCell.after_push(part, grid) = 
    ParticleInCell.discard!(part, grid) 

  Random.seed!(0)
  if isdefined(Main, :PROBLEM)
    problem = joinpath(pwd(), Main.PROBLEM)
  elseif length(ARGS) > 0
    problem = joinpath(pwd(), ARGS[1])
  else
    println("[!] No Main.PROBLEM nor ARGS defined")
    print("Provide problem filepath: ")
    filepath = readline()
    problem = joinpath(pwd(), filepath)
  end

  println("[*] Including problem definition from ", problem)
  include(problem)
  progress = Progress(ts)

  function run!()
    config = Config()
    config.grid = grid
    config.solver = solver
    config.pusher = pusher
    config.species = species
    config.interactions = interactions
    println("[*] Running simulation")
    @time ParticleInCell.solve(config, Δt, ts)
  end

  ParticleInCell.enter_loop() = begin
    start(Δt)
  end

  ParticleInCell.after_loop(i, t, dt) = begin
    update!(progress, i; showvalues = iteration(i, t, dt))
  end
    
  ParticleInCell.exit_loop() = stop()
end

julia_main() = IskraRuntime.run!()
julia_main()