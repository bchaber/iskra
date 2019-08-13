# iskra
A particle-in-cell simulator for spark simulations ("iskra" is the Polish word for spark).

# code
The code started as a loosely based version of Lubos Brieda's MATLAB code ES-PIC (https://www.particleincell.com/2011/particle-in-cell-example/). Right now only 2D (XY) electrostatic simulations of collisionless plasma is possible, but it will change in the future.

# running
To run this code few dependencies have to be installed in Julia (version >= 1.1). After that one has to add to `~/.julia/config/startup.jl` following lines:

```push!(LOAD_PATH, <PATH TO THE CLONED REPOSITORY>)```

After that change directory to `<PATH TO THE CLONED REPOSITORY>/problem` and run `julia problem.jl`.

# disclimer
Expect many sharp edges of this code. If you try to modify too much you will find a bear trying to eat you. Be careful.

