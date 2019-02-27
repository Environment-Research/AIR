
@defcomp climatedynamics begin

    ###############################################################################################################################
    # Model Parameters
    ###############################################################################################################################

    FORC    = Parameter(index=[time]) # Increase in radiative forcing (W/m²).
    fco22x  = Parameter()             # Forcings of equilibrium CO2 doubling (W/m²).
    t2xco2  = Parameter()             # Equilibrium temp impact (°C per doubling CO₂).
    tatm0   = Parameter()             # Initial atmospheric temperature anomaly (°C above pre-industrial).
    tocean0 = Parameter()             # Initial lower stratum temp change (°C above pre-industrial).
    c1      = Parameter()             # Climate equation coefficient for upper level.
    c3      = Parameter()             # Transfer coefficient upper to lower stratum.
    c4      = Parameter()             # Transfer coefficient for lower level.

    ###############################################################################################################################
    # Model Variables.
    ###############################################################################################################################
    TATM    = Variable(index=[time])  # Increase in temperature of atmosphere (°C above pre-industrial).
    TOCEAN  = Variable(index=[time])  # Increase in temperature of lower oceans (°C above pre-industrial).

    ###############################################################################################################################
    # Model Equations
    ###############################################################################################################################
    
    # Create run_timestep function to loop through each equation.
    function run_timestep(p, v, d, t)
        
        #Define function for atmospheric and ocean temperatures.
        if is_first(t)
            v.TATM[t] = p.tatm0
            v.TOCEAN[t] = p.tocean0
        else
            v.TATM[t] = v.TATM[t-1] + p.c1 * ((p.FORC[t] - (p.fco22x/p.t2xco2) * v.TATM[t-1]) - (p.c3 * (v.TATM[t-1] - v.TOCEAN[t-1])))
            v.TOCEAN[t] = v.TOCEAN[t-1] + p.c4 * (v.TATM[t-1] - v.TOCEAN[t-1])
        end
    end
end
