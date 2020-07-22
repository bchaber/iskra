# iskra
A particle-in-cell simulator for spark simulations ("iskra" is the Polish word for spark).

# origin
The code started as a loosely based version of Lubos Brieda's MATLAB code ES-PIC (https://www.particleincell.com/2011/particle-in-cell-example/). Right now it supports 1D (X) and 2D (XY) electrostatic simulations of plasma (2D simulations in cylindrical coordinate system will apear soon).
Another source of inspiration also comes from Lubos Brieda (https://github.com/particleincell/Starfish). It helped implement DSMC and MCC in iskra.

# features
In this code we try to keep a few separate modules: for the field solver, the grid, particle-in-cell algorithm and circuit simulations.

The circuit is coupled using surface charge density on the electrode' surface. For now, only RLC circuit can be connected between one of the electrode and the ground of the system. Charges hitting electrodes can be absorbed or reflected. The absorbed charge is accumulated and can interact with the circuit elements.

Iskra tries to keep working with diagnostics (particles' positions and other parameters, field distributions, etc.) as simple as possible. The computational code makes different quantities of interest available for saving them by the user. Currently, the diagnostics data are stored according to OpenPMD standard (http://openpmd.org) in HDF5 files. They **cannot** be interpreted out of the box using openpmd-viewer (due to the way HDF5.jl stores strings). You have to remove `.decode()` here and there from openpmd-viewer Python files...

# running
To run this code few dependencies have to be installed in Julia (version >= 1.1):

In Julia REPL press `]` and write the following commands to install iskra's depencies:
```
add DataStructures
add Unitful
add WriteVTK
```

After that one has to add to `~/.julia/config/startup.jl` following lines:

```push!(LOAD_PATH, <PATH TO THE CLONED REPOSITORY>)```

After that change directory to `<PATH TO THE CLONED REPOSITORY>/problem` and run `julia 01_single_electron.jl`.

# verification

This code has been partialy verified against XOOPIC from PTSG (http://ptsg.egr.msu.edu). We have compared two-stream instability simulation and got similar results from both codes. In the near future, DSMC and MCC modules will be verified.

# disclimer
Expect many sharp edges of this code. If you try to modify too much you will find a bear trying to eat you. Be careful.

