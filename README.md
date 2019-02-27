# RICE + Aerosol Impacts and Responses Model (RICE+AIR)

This repository contains source code for the RICE+AIR model presented in *The Impact of Human Health Co-benefits on Evaluations of Global Climate Policy, (Scovronick et al., Nature Communications, 2019)*.  RICE+AIR is written in the [Julia](www.julialang.org) programming language and uses the [Mimi framework](https://github.com/anthofflab/Mimi.jl) for integrated assessment model development.



## Software Requirements

RICE+AIR works with [Julia 1.0 and later versions](https://julialang.org/downloads/).

After downloading Julia, the required packages can be installed by executing the following commands using Julia's builtin package manager [Pkg](https://docs.julialang.org/en/v1/stdlib/Pkg/index.html):
````julia
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

This git repository uses a git submodule for the [Mimi implementation of RICE2010](https://github.com/anthofflab/mimi-rice-2010.jl). To ensure the submodule gets properly downloaded, make sure to use the git ``--recurse-submodules`` option when cloning the repository. If you cloned the repository without that option, you can issue the following two git commands to make sure the submodule is present on your system: ``git submodule init``, followed by ``git submodule update``.

## Issues?
If you have any trouble downloading the code or running RICE+AIR, please contact Frank Errickson at FrankErrickson@berkeley.edu
