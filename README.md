# RICE + Aerosol Impacts and Responses Model (RICE+AIR)

This repository contains source code for the RICE+AIR model presented in *The Impact of Human Health Co-benefits on Evaluations of Global Climate Policy, (Scovronick et al., Nature Communications, 2019)*.  RICE+AIR is written in the [Julia](www.julialang.org) programming language and uses the [Mimi framework](https://github.com/anthofflab/Mimi.jl) for integrated assessment model development.



## Software Requirements

RICE+AIR works with [Julia 1.1 and later versions](https://julialang.org/downloads/).  You will also need to make a one time connection to the [Mimi project's online model registry](https://github.com/mimiframework/MimiRegistry) (instructions below).

After downloading Julia, the required packages can be installed by executing the following commands using Julia's builtin package manager [Pkg](https://docs.julialang.org/en/v1/stdlib/Pkg/index.html):
````julia
using Pkg
Pkg.add("Mimi")
Pkg.add("NLopt")
Pkg.add("ExcelReaders")
Pkg.add("ExcelFiles")
Pkg.add("DataFrames")
Pkg.add("CSVFiles")
````


## Running RICE+AIR

(1) Open up Julia and set your current working directory to the **rice_air** folder.
````julia
cd("local_path_to_folder/rice_air")
````
(2) The **user_interface.jl** file in the **src** folder contains options that make it possible to run RICE+AIR with different model settings.  Users can change any value/setting under the *RICE+AIR Parameters to Change* and *Choices About Your Analysis & Optimzation* headings.  After making any changes, resave the **user_interface.jl** file.

(3) Next, run the **user_interface.jl** file.  This will evaluate RICE+AIR to find the optimal climate policy and will also save several key results in the **results** folder.
````julia
include("src/user_interface.jl")
````


## Note on Cloning the repository

This git repository uses the [Mimi implementation of RICE2010](https://github.com/anthofflab/MimiRICE2010.jl). The Mimi project maintains an [online registry](https://github.com/mimiframework/MimiRegistry) of Mimi models that includes RICE2010.  Once you connect your Julia installation to the Mimi registry, the models can be used in a similar manner to Julia packages (no longer requiring the use of git submodules for multi-model projects).

To connect to the Mimi registry, you need to run the following command only once using the Julia package REPL (note, you can enter the package REPL using the `]` key and exit using the `backspace` key):
````julia
registry add https://github.com/mimiframework/MimiRegistry.git
````

## Issues?
If you have any trouble downloading the code or running RICE+AIR, please contact Frank Errickson at FrankErrickson@berkeley.edu
