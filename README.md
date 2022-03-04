<p align="center">
  <img src="https://github.com/bchaber/iskra/blob/master/logo.svg" width="256px" alt="The logo is based on a symbol of the highest Slavic god: Perun found at https://commons.wikimedia.org/wiki/File:Thundermarks.svg"/>
  <h1>iskra <a href="https://zenodo.org/badge/latestdoi/202223612"><img src="https://zenodo.org/badge/202223612.svg" alt="DOI"></a></h1>
</p>

A particle-in-cell simulator for spark simulations ("iskra" is the Polish word for spark).

![](../../blob/master/img/two-stream-rho-evolution.png)

# :construction_worker: announcement
Currently the code is being rewritten and will be significantly faster. The field solvers were moved to a separate repositories (namely: [pole](https://github.com/bchaber/pole) and [fala](https://github.com/bchaber/fala)). We plan to put the improved version of the code in Spring 2022.

# origin
The code started as a loosely based version of Lubos Brieda's MATLAB code ES-PIC (https://www.particleincell.com/2011/particle-in-cell-example/). Right now it supports 1D (X) and 2D (XY) electrostatic simulations of plasma (2D simulations in cylindrical coordinate system will appear soon).
Another source of inspiration also comes from Lubos Brieda (https://github.com/particleincell/Starfish). It helped implement DSMC and MCC in iskra.

# publications
Below is the list of publications using this project:

[1] B. Chaber, "Particle-in-Cell code for gas discharge simulations," 2020 IEEE 21st International Conference on Computational Problems of Electrical Engineering (CPEE), 2020, pp. 1-4, doi: 10.1109/CPEE50798.2020.9238682.
[2] W. Łodyga, B. Chaber, "Parallel implementation of a Particle-in-Cell code in Julia programming language," 2021 22nd International Conference on Computational Problems of Electrical Engineering (CPEE), 2021, doi: 10.1109/CPEE54040.2021.9585274
[3] B. Chaber, "Cell-centered, geometric multigrid field solver for Particle-in-Cell simulations," 2021 22nd International Conference on Computational Problems of Electrical Engineering (CPEE), 2021, doi: 10.1109/CPEE54040.2021.9585258

# features
Currently, we are focused on the basic Particle-in-Cell implementation for a collissionless plasma.
We have dropped OpenPMD support for now, as we can't make it fast enough at the moment.

# running
To run the project install the IJulia:

```
$ cd iskra/
$ ls *.toml
Manifest.toml	Project.toml
$ julia --project
julia> import Pkg
julia> Pkg.add("IJulia")
julia> import IJulia
julia> IJulia.notebook()
```

Alternatively, you start Jupyter Notebook with: `$ julia --project -e "import IJulia; IJulia.notebook()"`

# verification

This code has been partialy verified against XOOPIC from PTSG (http://ptsg.egr.msu.edu).
We have compared two-stream instability simulation and got similar results from both codes.

# disclimer
Expect many sharp edges of this code. If you try to modify too much you will find a bear trying to eat you. Be careful.

