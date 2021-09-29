<p align="center">
  <img src="https://github.com/bchaber/iskra/blob/master/logo.svg" width="256px" alt="The logo is based on a symbol of the highest Slavic god: Perun found at https://commons.wikimedia.org/wiki/File:Thundermarks.svg"/>
  <h1>iskra</h1>
</p>

A particle-in-cell simulator for spark simulations ("iskra" is the Polish word for spark).

![](../../blob/master/img/two-stream-rho-evolution.png)

# origin
The code started as a loosely based version of Lubos Brieda's MATLAB code ES-PIC (https://www.particleincell.com/2011/particle-in-cell-example/). Right now it supports 1D (X) and 2D (XY) electrostatic simulations of plasma (2D simulations in cylindrical coordinate system will appear soon).
Another source of inspiration also comes from Lubos Brieda (https://github.com/particleincell/Starfish). It helped implement DSMC and MCC in iskra.

# publications
Below is the list of publications using this project:

[1] B. Chaber, "Particle-in-Cell code for gas discharge simulations," 2020 IEEE 21st International Conference on Computational Problems of Electrical Engineering (CPEE), 2020, pp. 1-4, doi: 10.1109/CPEE50798.2020.9238682.

# features
In this code we try to keep a few separate modules: for the field solver, the grid, particle-in-cell algorithm and circuit simulations.

The circuit is coupled using surface charge density on the electrode' surface. For now, only RLC circuit can be connected between one of the electrode and the ground of the system. Charges hitting electrodes can be absorbed or reflected. The absorbed charge is accumulated and can interact with the circuit elements.

Iskra tries to keep working with diagnostics (particles' positions and other parameters, field distributions, etc.) as simple as possible. The computational code makes different quantities of interest available for saving them by the user. Currently, the diagnostics data are stored according to OpenPMD standard (http://openpmd.org) in HDF5 files. They **cannot** be interpreted out of the box using openpmd-viewer (due to the way HDF5.jl stores strings). You have to remove `.decode()` here and there from openpmd-viewer Python files...

# running
To run this code few dependencies have to be installed in Julia (version >= 1.5):

You can use the environment of the package:

```
$ cd iskra/
$ ls *.toml
Manifest.toml	Project.toml
$ julia --project src/iskra.jl problem/01_single_electron.jl
...
```

Alternatively, you can run it from REPL:

```
$ julia --project
julia> PROBLEM = "problem/01_single_electron.jl"
julia> include("src/iskra.jl")
...
```

Please, mind that the during the saving process of diagnostics data the current
directory is changed, so to re-run some simulation you might want to change the
directory in REPL.

# verification

This code has been partialy verified against XOOPIC from PTSG (http://ptsg.egr.msu.edu). We have compared two-stream instability simulation and got similar results from both codes. In the near future, DSMC and MCC modules will be verified.

![](../../blob/master/img/rf-discharge-density-evolution.png)

# disclimer
Expect many sharp edges of this code. If you try to modify too much you will find a bear trying to eat you. Be careful.

