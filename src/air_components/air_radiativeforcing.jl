
@defcomp air_radiativeforcing begin

    ###############################################################################################################################
    # Model Indices
    ###############################################################################################################################
    regions   = Index()                                            # 12 regions in RICE2010.

    ###############################################################################################################################
    # Model Parameters
    ###############################################################################################################################
    fco22x                      = Parameter()                      # Forcings of equilibrium CO2 doubling (W/m²).
    mat1                        = Parameter()                      # 2015 atmospheric CO₂ concentration (from RICE2010).
    MAT                         = Parameter(index=[time])          # Carbon concentration increase in atmosphere (GtC from 1750)
    MATSUM                      = Parameter(index=[time])          # Sum of MAT[t] and MAT[t+1] to use in FORC[t] for radiativeforcing component (follows RICE2010).
    forcoth                     = Parameter(index=[time])          # Exogenous forcing for other greenhouse gases
    u0_SO2                      = Parameter(index=[regions])       # Additive term for calculating SO₂ r-coefficient.
    u0_NOX                      = Parameter(index=[regions])       # Additive term for calculating NOx r-coefficient.
    u0_BC                       = Parameter(index=[regions])       # Additive term for calculating BC r-coefficient.
    u0_OC                       = Parameter(index=[regions])       # Additive term for calculating OC r-coefficient.
    u1_SO2                      = Parameter(index=[regions])       # Multiplicative term for calculating SO₂ r-coefficient.
    u1_NOX                      = Parameter(index=[regions])       # Multiplicative term for calculating NOx r-coefficient.
    u1_BC                       = Parameter(index=[regions])       # Multiplicative term for calculating BC r-coefficient.
    u1_OC                       = Parameter(index=[regions])       # Multiplicative term for calculating OC r-coefficient.
    base_emissions_SO2          = Parameter(index=[time, regions]) # Baseline (pre-mitigation) SO₂ emissions.
    base_emissions_NOX          = Parameter(index=[time, regions]) # Baseline (pre-mitigation) NOx emissions.
    base_emissions_BC           = Parameter(index=[time, regions]) # Baseline (pre-mitigation) BC emissions.
    base_emissions_OC           = Parameter(index=[time, regions]) # Baseline (pre-mitigation) OC emissions.
    Δemissions_SO2              = Parameter(index=[time, regions]) # Change (co-reduction) in SO₂ emissions due to CO₂ mitigation.
    Δemissions_NOX              = Parameter(index=[time, regions]) # Change (co-reduction) in NOx emissions due to CO₂ mitigation.
    Δemissions_BC               = Parameter(index=[time, regions]) # Change (co-reduction) in BC emissions due to CO₂ mitigation.
    Δemissions_OC               = Parameter(index=[time, regions]) # Change (co-reduction) in OC emissions due to CO₂ mitigation.

    ###############################################################################################################################
    # Model Variables.
    ###############################################################################################################################
    base_aerosol_rf             = Variable(index=[time])           # Radiative forcing from baseline aerosol emission levels (W/m²).
    Δaerosol_rf                 = Variable(index=[time])           # Change in raditive forcing due to aerosol co-reduction (W/m²).
    post_mitigation_aerosol_rf  = Variable(index=[time])           # Total radiative forcing due to aerosols after accounting for aerosol co-reductions (W\m²).
    FORC                        = Variable(index=[time])           # Increase in radiative forcing (RICE2010 has units of W/m² from 1900).
    r_coefficient_SO2           = Variable(index=[time, regions])  # Coefficient to relate regional change in SO₂ emissions to change in global radiative forcing.
    r_coefficient_NOX           = Variable(index=[time, regions])  # Coefficient to relate regional change in NOx emissions to change in global radiative forcing.
    r_coefficient_BC            = Variable(index=[time, regions])  # Coefficient to relate regional change in BC emissions to change in global radiative forcing.
    r_coefficient_OC            = Variable(index=[time, regions])  # Coefficient to relate regional change in OC emissions to change in global radiative forcing.
    base_rf_SO2                 = Variable(index=[time, regions])  # Radiative forcing from baseline SO₂ emission levels (W/m²).
    base_rf_NOX                 = Variable(index=[time, regions])  # Radiative forcing from baseline NOx emission levels (W/m²).
    base_rf_BC                  = Variable(index=[time, regions])  # Radiative forcing from baseline BC emission levels (W/m²).
    base_rf_OC                  = Variable(index=[time, regions])  # Radiative forcing from baseline OC emission levels (W/m²).
    Δrf_SO2                     = Variable(index=[time, regions])  # Change in raditive forcing due to SO₂ co-reduction (W/m²).
    Δrf_NOX                     = Variable(index=[time, regions])  # Change in raditive forcing due to NOx co-reduction (W/m²).
    Δrf_BC                      = Variable(index=[time, regions])  # Change in raditive forcing due to BC co-reduction (W/m²).
    Δrf_OC                      = Variable(index=[time, regions])  # Change in raditive forcing due to OC co-reduction (W/m²).
    co2_forcing                 = Variable(index=[time])           # Radiative forcing from CO₂ (W/m²).

    ###############################################################################################################################
    # Model Equations
    ###############################################################################################################################
    
    # Create run_timestep function to loop through each equation.
    function run_timestep(p, v, d, t)

        # Time dependency term used to calculate r-coefficient.
        time_term = 1995.0 + t.t * 10.0

        for r in d.regions

            # Calculate r-coefficient for each aerosol type over time.
            v.r_coefficient_SO2[t,r] = p.u1_SO2[r] * time_term + p.u0_SO2[r]
            v.r_coefficient_NOX[t,r] = p.u1_NOX[r] * time_term + p.u0_NOX[r]
            v.r_coefficient_BC[t,r]  = p.u1_BC[r]  * time_term + p.u0_BC[r]
            v.r_coefficient_OC[t,r]  = p.u1_OC[r]  * time_term + p.u0_OC[r]

            # Calculate baseline radiative forcing for each aerosol type.
            v.base_rf_SO2[t,r] = v.r_coefficient_SO2[t,r] * p.base_emissions_SO2[t,r]
            v.base_rf_NOX[t,r] = v.r_coefficient_NOX[t,r] * p.base_emissions_NOX[t,r]
            v.base_rf_BC[t,r]  = v.r_coefficient_BC[t,r]  * p.base_emissions_BC[t,r]
            v.base_rf_OC[t,r]  = v.r_coefficient_OC[t,r]  * p.base_emissions_OC[t,r]

            # Calculate change in radiative forcing due to co-reduction for each aerosol type.
            v.Δrf_SO2[t,r] = v.r_coefficient_SO2[t,r] * p.Δemissions_SO2[t,r]
            v.Δrf_NOX[t,r] = v.r_coefficient_NOX[t,r] * p.Δemissions_NOX[t,r]
            v.Δrf_BC[t,r]  = v.r_coefficient_BC[t,r]  * p.Δemissions_BC[t,r]
            v.Δrf_OC[t,r]  = v.r_coefficient_OC[t,r]  * p.Δemissions_OC[t,r]
        end

        # Calculate total baseline, change in, and post mitigation radiative forcing for all aerosol types.
        v.base_aerosol_rf[t] = sum(v.base_rf_SO2[t,:] .+ v.base_rf_NOX[t,:] .+ v.base_rf_BC[t,:] .+ v.base_rf_OC[t,:])
        v.Δaerosol_rf[t]     = sum(v.Δrf_SO2[t,:] .+ v.Δrf_NOX[t,:] .+ v.Δrf_BC[t,:] .+ v.Δrf_OC[t,:])
        v.post_mitigation_aerosol_rf[t] = v.base_aerosol_rf[t] - v.Δaerosol_rf[t]

        # Calculate global radiative forcing (accounting for CO₂, exogenous RF sources, and endogenous aerosol effect).
        if is_first(t)
            v.co2_forcing[t] = p.fco22x * ((log10((((p.MAT[t] + p.mat1)/2)+0.000001)/596.4)/log10(2)))
            v.FORC[t] = v.co2_forcing[t] + p.forcoth[t] + v.post_mitigation_aerosol_rf[t]
        else
            v.co2_forcing[t] = p.fco22x * ((log10((((p.MATSUM[t])/2)+0.000001)/596.4)/log10(2)))
            v.FORC[t] =  v.co2_forcing[t] + p.forcoth[t] + v.post_mitigation_aerosol_rf[t]
        end
    end
end
