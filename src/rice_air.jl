
# Load files needed to construct RICE+AIR.
include("air_components/air_baseline_emissions.jl")
include("air_components/air_coreduction.jl")
include("air_components/air_lifeyears.jl")
include("air_components/air_radiativeforcing.jl")
include("air_components/air_consumption.jl")
include("modified_rice_components/welfare_modified.jl")
include("rice2010/src/rice2010.jl")
include("modified_rice_components/co2cycle_modified.jl")
include("modified_rice_components/climatedynamics_modified.jl")
include("helper_functions.jl")

# Load Rice2010 module.
using .Rice2010


####################################################################################################
# Create a function to construct RICE+AIR given user defined parameter settings.
####################################################################################################

function construct_rice_air(inputs::RICE_AIR_inputs)

    # Read in RICE+AIR Parameters and data.
    riceair_params = load_riceair_params(inputs.SSP_scenario)

    # Read in RICE2010 parameters.
    rice_params = Rice2010.getrice2010parameters("data/RICE_2010_base_000.xlsm")

    # Set RICE2010 population for initial value, then use updated UN population values and rescale to correct units.
    rice_air_population = zeros(inputs.nsteps,12)
    rice_air_population[1,:] = rice_params[:l][1,:]
    rice_air_population[2:end,:] = (riceair_params[:population][1:(inputs.nsteps-1),:])./1000.0

    # Create an instance of RICE2010, and then add in AIR components.
    m = Rice2010.getrice()

    # Set index for 5 aerosol types in RICE+AIR.
    set_dimension!(m, :aerosols, ["SO2", "PM2.5", "NOX", "BC", "OC"])

    # Delete RICE2010 Radiative Forcing, CO₂ Cycle, Climate, and Welfare components.
    delete!(m, :co2cycle)
    delete!(m, :climatedynamics)
    delete!(m, :radiativeforcing)
    delete!(m, :welfare)

    # Add RICE+AIR specific components to remaining RICE2010 model structure.
    add_comp!(m, co2cycle,               after=:emissions)
    add_comp!(m, air_baseline_emissions, after=:co2cycle)
    add_comp!(m, air_coreduction,        after=:air_baseline_emissions)
    add_comp!(m, air_radiativeforcing,   after=:air_coreduction)
    add_comp!(m, climatedynamics,        after=:air_radiativeforcing)
    add_comp!(m, air_lifeyears,          after=:neteconomy)
    add_comp!(m, air_consumption,        after=:air_lifeyears)
    add_comp!(m, welfare,                after=:air_consumption)

    # -------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------
    # Set all model parameters for new RICE+AIR components.
    # -------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------

    #---- GROSS ECONOMIC OUTPUT ----#
    set_param!(m, :grosseconomy, :l, rice_air_population)

    #---- NET ECONOMIC OUTPUT ----#
    set_param!(m, :neteconomy, :l, rice_air_population)

    #---- EMISSIONS ----#
    set_param!(m, :emissions, :MIU, zeros(inputs.nsteps, 12))

    #---- CO₂ CYCLE ----#
    set_param!(m, :co2cycle, :mat0, rice_params[:mat0])
    set_param!(m, :co2cycle, :mu0, rice_params[:mu0])
    set_param!(m, :co2cycle, :ml0, rice_params[:ml0])
    set_param!(m, :co2cycle, :b12, rice_params[:b12])
    set_param!(m, :co2cycle, :b23, rice_params[:b23])
    set_param!(m, :co2cycle, :b11, rice_params[:b11])
    set_param!(m, :co2cycle, :b21, rice_params[:b21])
    set_param!(m, :co2cycle, :b22, rice_params[:b22])
    set_param!(m, :co2cycle, :b32, rice_params[:b32])
    set_param!(m, :co2cycle, :b33, rice_params[:b33])

    #---- BASELINE AEROSOL EMISSIONS ----#
    set_param!(m, :air_baseline_emissions, :population, rice_air_population)
    set_param!(m, :air_baseline_emissions, :SO2₀, riceair_params[:SO2₀])
    set_param!(m, :air_baseline_emissions, :PM25₀, riceair_params[:PM25₀])
    set_param!(m, :air_baseline_emissions, :NOX₀, riceair_params[:NOX₀])
    set_param!(m, :air_baseline_emissions, :omega_1, riceair_params[:omega_1])
    set_param!(m, :air_baseline_emissions, :omega_2, riceair_params[:omega_2])
    set_param!(m, :air_baseline_emissions, :phi_1, riceair_params[:phi_1])
    set_param!(m, :air_baseline_emissions, :phi_2, riceair_params[:phi_2])
    set_param!(m, :air_baseline_emissions, :phi_3, riceair_params[:phi_3])
    set_param!(m, :air_baseline_emissions, :K, inputs.kuznets_term)

    #---- AEROSOL CO-REDUCTIONS FROM SHARED SOURCES ----#
    set_param!(m, :air_coreduction, :MIU, zeros(inputs.nsteps, 12))
    set_param!(m, :air_coreduction, :κ, riceair_params[:kappa])

    #---- RADIATIVE FORCING WITH ENDOGENOUS AEROSOLS ----#
    set_param!(m, :air_radiativeforcing, :u0_SO2, riceair_params[:u0_SO2])
    set_param!(m, :air_radiativeforcing, :u0_NOX, riceair_params[:u0_NOX])
    set_param!(m, :air_radiativeforcing, :u0_BC,  riceair_params[:u0_BC])
    set_param!(m, :air_radiativeforcing, :u0_OC,  riceair_params[:u0_OC])
    set_param!(m, :air_radiativeforcing, :u1_SO2, riceair_params[:u1_SO2])
    set_param!(m, :air_radiativeforcing, :u1_NOX, riceair_params[:u1_NOX])
    set_param!(m, :air_radiativeforcing, :u1_BC,  riceair_params[:u1_BC])
    set_param!(m, :air_radiativeforcing, :u1_OC,  riceair_params[:u1_OC])
    set_param!(m, :air_radiativeforcing, :forcoth, riceair_params[:total_exogeous_rf])
    set_param!(m, :air_radiativeforcing, :fco22x, rice_params[:fco22x])
    set_param!(m, :air_radiativeforcing, :mat1, rice_params[:mat1])

    #---- CLIMATE MODEL ----#
    set_param!(m, :climatedynamics, :fco22x, rice_params[:fco22x])
    set_param!(m, :climatedynamics, :t2xco2, rice_params[:t2xco2])
    set_param!(m, :climatedynamics, :tatm0, rice_params[:tatm0])
    set_param!(m, :climatedynamics, :tocean0, rice_params[:tocean0])
    set_param!(m, :climatedynamics, :c1, rice_params[:c1])
    set_param!(m, :climatedynamics, :c3, rice_params[:c3])
    set_param!(m, :climatedynamics, :c4, rice_params[:c4])
 
    #---- LIFEYEARS GAINED FROM CO-REDUCTIONS ----#
    set_param!(m, :air_lifeyears, :τ, inputs.tau)
    set_param!(m, :air_lifeyears, :upper_threshold, 10000.0)
    set_param!(m, :air_lifeyears, :β, riceair_params[:β])
    set_param!(m, :air_lifeyears, :exposure₀, riceair_params[:exposure_2005])
    set_param!(m, :air_lifeyears, :source_receptor_SO2, riceair_params[:sr_SO2])
    set_param!(m, :air_lifeyears, :source_receptor_PM25, riceair_params[:sr_PM25])
    set_param!(m, :air_lifeyears, :source_receptor_NOX, riceair_params[:sr_NOX])
    set_param!(m, :air_lifeyears, :Θ , riceair_params[:Θ])
    set_param!(m, :air_lifeyears, :SO2₀, riceair_params[:SO2₀])
    set_param!(m, :air_lifeyears, :PM25₀, riceair_params[:PM25₀])
    set_param!(m, :air_lifeyears, :NOX₀, riceair_params[:NOX₀])
    set_param!(m, :air_lifeyears, :base_deaths , riceair_params[:base_deaths])

    #---- MONTEIZED HEALTH BENEFITS & CONSUMPTION ----#
    set_param!(m, :air_consumption, :Hyears, inputs.Hyears)
    set_param!(m, :air_consumption, :VOLY_elasticity, inputs.VOLY_elasticity)
    set_param!(m, :air_consumption, :pop, rice_air_population)
    set_param!(m, :air_consumption, :use_VSL, inputs.use_VSL)

    #---- ECONOMIC WELFARE ----#
    set_param!(m, :welfare, :pop, rice_air_population)
    set_param!(m, :welfare, :rho, inputs.rho)
    set_param!(m, :welfare, :eta, inputs.eta)

    # -------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------
    # Create connections between coupled RICE+AIR components.
    # -------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------

    #---- BASELINE AEROSOL EMISSIONS ----#
    connect_param!(m, :air_baseline_emissions, :gross_output, :grosseconomy, :YGROSS)

    #---- AEROSOL CO-REDUCTIONS FROM SHARED SOURCES ----#
    connect_param!(m, :air_coreduction, :base_emissions_SO2,  :air_baseline_emissions, :base_emissions_SO2)
    connect_param!(m, :air_coreduction, :base_emissions_PM25, :air_baseline_emissions, :base_emissions_PM25)
    connect_param!(m, :air_coreduction, :base_emissions_NOX,  :air_baseline_emissions, :base_emissions_NOX)
    connect_param!(m, :air_coreduction, :base_emissions_BC,   :air_baseline_emissions, :base_emissions_BC)
    connect_param!(m, :air_coreduction, :base_emissions_OC,   :air_baseline_emissions, :base_emissions_OC)

    #---- CO₂ CYCLE ----#
    connect_param!(m, :co2cycle, :E, :emissions, :E)

    #---- LIFEYEARS GAINED FROM CO-REDUCTIONS ----#
    connect_param!(m, :air_lifeyears, :base_emissions_SO2,  :air_baseline_emissions, :base_emissions_SO2)
    connect_param!(m, :air_lifeyears, :base_emissions_PM25, :air_baseline_emissions, :base_emissions_PM25)
    connect_param!(m, :air_lifeyears, :base_emissions_NOX,  :air_baseline_emissions, :base_emissions_NOX)
    connect_param!(m, :air_lifeyears, :Δemissions_SO2,      :air_coreduction,        :Δemissions_SO2)
    connect_param!(m, :air_lifeyears, :Δemissions_PM25,     :air_coreduction,        :Δemissions_PM25)
    connect_param!(m, :air_lifeyears, :Δemissions_NOX,      :air_coreduction,        :Δemissions_NOX)

    #---- RADIATIVE FORCING WITH ENDOGENOUS AEROSOLS ----#
    connect_param!(m, :air_radiativeforcing, :MAT, :co2cycle, :MAT)
    connect_param!(m, :air_radiativeforcing, :MATSUM, :co2cycle, :MATSUM)
    connect_param!(m, :air_radiativeforcing, :base_emissions_SO2,  :air_baseline_emissions, :base_emissions_SO2)
    connect_param!(m, :air_radiativeforcing, :base_emissions_NOX,  :air_baseline_emissions, :base_emissions_NOX)
    connect_param!(m, :air_radiativeforcing, :base_emissions_BC,  :air_baseline_emissions, :base_emissions_BC)
    connect_param!(m, :air_radiativeforcing, :base_emissions_OC,  :air_baseline_emissions, :base_emissions_OC)
    connect_param!(m, :air_radiativeforcing, :Δemissions_SO2,      :air_coreduction,        :Δemissions_SO2)
    connect_param!(m, :air_radiativeforcing, :Δemissions_NOX,      :air_coreduction,        :Δemissions_NOX)
    connect_param!(m, :air_radiativeforcing, :Δemissions_BC,      :air_coreduction,        :Δemissions_BC)
    connect_param!(m, :air_radiativeforcing, :Δemissions_OC,      :air_coreduction,        :Δemissions_OC)

    #---- CLIMATE MODEL ----#
    connect_param!(m, :climatedynamics, :FORC, :air_radiativeforcing, :FORC)

    #---- SEA LEVEL RISE ----#
    connect_param!(m, :sealevelrise, :TATM, :climatedynamics, :TATM)

    #---- ECONOMIC CLIMATE DAMAGES ----#
    connect_param!(m, :damages, :TATM, :climatedynamics, :TATM)

    #---- MONTEIZED HEALTH BENEFITS & CONSUMPTION ----#
    connect_param!(m, :air_consumption, :CPC, :neteconomy, :CPC)
    connect_param!(m, :air_consumption, :lifeyears, :air_lifeyears, :lifeyears)
    connect_param!(m, :air_consumption, :avoided_deaths, :air_lifeyears, :avoided_deaths)
    connect_param!(m, :air_consumption, :Y, :neteconomy, :Y)

    #---- ECONOMIC WELFARE ----#
    connect_param!(m, :welfare, :CPC_with_health, :air_consumption, :CPC_with_health)

    # Return RICE+AIR model.
    return m
end
