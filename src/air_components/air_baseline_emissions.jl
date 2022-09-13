
@defcomp air_baseline_emissions begin

    ###############################################################################################################################
    # Model Indices
    ###############################################################################################################################
    regions  = Index()                                              # 12 regions in RICE2010.
    aerosols = Index()                                              # 5 aerosols (SO₂, PM2.5, NOx, BC, OC).

    ###############################################################################################################################
    # Model Parameters
    ###############################################################################################################################
    K                   = Parameter()                               # Kuznet's scaling term for rate of clean up.
    SO2₀                = Parameter(index=[regions])                # Initial value for baseline (pre-mitigation) SO₂ emissions (Gg).
    PM25₀               = Parameter(index=[regions])                # Initial value for baseline (pre-mitigation) PM2.5 emissions (Gg).
    NOX₀                = Parameter(index=[regions])                # Initial value for baseline (pre-mitigation) NOx emissions (Gg).
    population          = Parameter(index=[time, regions])          # Regional population levels (millions of people).
    gross_output        = Parameter(index=[time, regions])          # Gross economic output ($trillions).
    omega_1             = Parameter(index=[regions, aerosols])      # Fitting parameter for emissions intensity factor equation.
    omega_2             = Parameter(index=[regions, aerosols])      # Fitting parameter for emissions intensity factor equation.
    phi_1               = Parameter(index=[regions, aerosols])      # Fitting parameter for emissions intensity factor equation.
    phi_2               = Parameter(index=[regions, aerosols])      # Fitting parameter for emissions intensity factor equation.
    phi_3               = Parameter(index=[regions, aerosols])      # Fitting parameter for emissions intensity factor equation.
    
    ###############################################################################################################################
    # Model Variables
    ###############################################################################################################################
    kuznets_income      = Variable(index=[time, regions])           # Kuznet's income scaling term to adjust autonomous rate of air pollution clean up.
    base_emissions_SO2  = Variable(index=[time, regions])           # Baseline (pre-mitigation) SO₂ emissions (Gg).
    base_emissions_PM25 = Variable(index=[time, regions])           # Baseline (pre-mitigation) PM2.5 emissions (Gg).
    base_emissions_NOX  = Variable(index=[time, regions])           # Baseline (pre-mitigation) NOx emissions (Gg).
    base_emissions_BC   = Variable(index=[time, regions])           # Baseline (pre-mitigation) BC emissions (Gg).
    base_emissions_OC   = Variable(index=[time, regions])           # Baseline (pre-mitigation) OC emissions (Gg).
    e_intensity         = Variable(index=[time, regions, aerosols]) # Emission intensity factors (ordering goes SO₂, PM2.5, NOx, BC, OC).

    ###############################################################################################################################
    # Model Equations
    ###############################################################################################################################
    
    # Create run_timestep function to loop through each equation.
    function run_timestep(p, v, d, t)

        for r in d.regions

            # Calculate Kuznets relevant income.
            v.kuznets_income[t,r] = (p.gross_output[TimestepIndex(1),r] / (p.population[TimestepIndex(1),r]/1000.0)) + p.K * ((p.gross_output[t,r] / (p.population[t,r]/1000.0)) - (p.gross_output[TimestepIndex(1),r] / (p.population[TimestepIndex(1),r]/1000.0)))

            for a in d.aerosols
                # Calculate emission intensity factors for the five aerosols.
                v.e_intensity[t,r,a] = p.phi_1[r,a] * (exp(-p.omega_1[r,a] * v.kuznets_income[t,r]) + p.phi_2[r,a] * exp(-p.omega_2[r,a] * v.kuznets_income[t,r])) + p.phi_3[r,a]
            end

            # Set/calculate initial values and baseline (pre-mitigation) emissions for all five aerosols.
            if is_first(t)
                v.base_emissions_SO2[t,r]  = p.SO2₀[r]
                v.base_emissions_PM25[t,r] = p.PM25₀[r]
                v.base_emissions_NOX[t,r]  = p.NOX₀[r]
                v.base_emissions_BC[t,r]   = v.e_intensity[t,r,4] * p.gross_output[t,r]
                v.base_emissions_OC[t,r]   = v.e_intensity[t,r,5] * p.gross_output[t,r]
            else
                v.base_emissions_SO2[t,r]  = v.e_intensity[t,r,1] * p.gross_output[t,r]
                v.base_emissions_PM25[t,r] = v.e_intensity[t,r,2] * p.gross_output[t,r]
                v.base_emissions_NOX[t,r]  = v.e_intensity[t,r,3] * p.gross_output[t,r]
                v.base_emissions_BC[t,r]   = v.e_intensity[t,r,4] * p.gross_output[t,r]
                v.base_emissions_OC[t,r]   = v.e_intensity[t,r,5] * p.gross_output[t,r]
            end
        end
    end
end
