# Load required Julia packages.
using NLopt
using MimiRICE2010
using Mimi
using ExcelReaders
using ExcelFiles
using DataFrames
using CSVFiles

# Load RICE+AIR source code.
include("helper_functions.jl")
include("rice_air.jl")

####################################################################################################
# RICE+AIR PARAMETERS TO CHANGE
####################################################################################################

# The number of time periods to run the model for (defaults to 60, representing 2005-2595 with 10 year time steps).
nsteps = 60

# Pure rate of time preference.
rho = 0.015

# Elasticity of marginal utility of consumption.
eta = 1.5

# Lower exposure threshold below which no health benefits accrue (μg/m³).
τ = 5.8

# Kuznet's income scaling term to adjust autonomous rate of air pollution clean up.
K = 1.0

# SSP Scenario to estimate effectiveness of CO₂ mitigation in co-reducing aerosol emissions (options = :SSP1 - :SSP5)
SSP_scenario = :SSP2

# Should the health benefit be calculated as avoided deaths with a VSL? (true = VSL valuation, false = lifeyears valuation).
use_VSL = false

# Number of years of per capita consumption to use in VOLY calculations for each region.
Hyears = 2.0

# Consumption elasticity term on the VOLY (defaults to 1.0).
VOLY_elasticity = 1.0

####################################################################################################
# CHOICES ABOUT YOUR ANALYSIS & OPTIMZATION
####################################################################################################

# Name of folder to store your results in (a folder will be created with this name).
results_folder = "My Results"

# Do you also want to perform a reference case optimization run (with health co-benefits turned off)?
run_reference_case = true

# Number of 10-year timesteps to find optimal carbon tax for (after which model assumes full decarbonization). 
n_objective = 25

# Optimization algorithm (:symbol). See options at http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
opt_algorithm = :LN_SBPLX

# Maximum time in seconds to run (in case optimization does not converge).
stop_time = 1000

# Relative tolerance criteria for convergence (will stop if |Δf| / |f| < tolerance from one iteration to the next.)
tolerance = 1e-14


####################################################################################################
####################################################################################################
# RUN EVERYTHING & SAVE KEY RESULTS
####################################################################################################
####################################################################################################

# Create a folder to save results.
output_directory = joinpath((@__DIR__), "../results", results_folder)
mkdir(output_directory)

# Select optimization starting point as a fraction of maximum carbon tax in each period.
rice_parameters = MimiRICE2010.getrice2010parameters("data/RICE_2010_base_000.xlsm")
starting_vals = maximum(rice_parameters[:pbacktime], dims=2)[2:(n_objective+1)] .* 1000 .* 0.001

# Combine user-defined inputs to create a particular specification of RICE+AIR.
inputs = RICE_AIR_inputs(nsteps, rho, eta, τ, K, SSP_scenario, Hyears, use_VSL, VOLY_elasticity)


#---------------------------------------------------------------------------------------------------
# Optimization Run for RICE+AIR with Health Co-Benefits
#---------------------------------------------------------------------------------------------------

# Optimize RICE+AIR (return regional decarbonization rates, optimal carbon tax, an instance of RICE+AIR, and optimization output vector).
cobenefit_decarbonization, cobenefit_tax, cobenefit_model, minx_cobenefit = optimize_rice_air(inputs, opt_algorithm, n_objective, stop_time, tolerance, rice_parameters[:pbacktime], true, starting_vals)

# Run instance of model with optimal decarbonization policy.
update_param!(cobenefit_model, :MIU, cobenefit_decarbonization)
run(cobenefit_model)

# Save key model output (global temperature anomaly, carbon tax, global decarbonization rate).
save(joinpath(output_directory, "cobenefit_temperature.csv"), DataFrame(temperature=cobenefit_model[:climatedynamics, :TATM]))
save(joinpath(output_directory, "cobenefit_carbon_tax.csv"), DataFrame(tax=cobenefit_tax))
save(joinpath(output_directory, "cobenefit_global_mitigation.csv"), DataFrame(mitigation=global_mitigation_mix(cobenefit_model)[:,1]))

#---------------------------------------------------------------------------------------------------
# Optimization Run for RICE+AIR Reference Case (no Health Co-Benefits)
#---------------------------------------------------------------------------------------------------

if run_reference_case

    # Optimize Reference Case (return regional decarbonization rates, optimal carbon tax, an instance of RICE+AIR, and optimization output vector).
    reference_decarbonization, reference_tax, reference_model, minx_reference = optimize_rice_air(inputs, opt_algorithm, n_objective, stop_time, tolerance, rice_parameters[:pbacktime], false, starting_vals)

    # Run instance of model with optimal decarbonization policy for reference case.
    update_param!(reference_model, :MIU, reference_decarbonization)
    run(reference_model)

    # Save key model output for reference case (global temperature anomaly, carbon tax, global decarbonization rate).
    save(joinpath(output_directory, "reference_temperature.csv"), DataFrame(temperature=reference_model[:climatedynamics, :TATM]))
    save(joinpath(output_directory, "reference_carbon_tax.csv"), DataFrame(tax=reference_tax))
    save(joinpath(output_directory, "reference_global_mitigation.csv"), DataFrame(mitigation=global_mitigation_mix(reference_model)[:,1]))

end

println("Finished with all RICE+AIR optimizations.")
