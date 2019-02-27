
@defcomp air_lifeyears begin

    ###############################################################################################################################
    # Model Indices
    ###############################################################################################################################
    regions  = Index()                                       # 12 regions in RICE2010.
    aerosols = Index()                                       # 5 aerosols (SO₂, PM2.5, NOx, BC, OC).
    
    ###############################################################################################################################
    # Model Parameters
    ###############################################################################################################################
    τ                     = Parameter()                      # Exposure threshold below which no health benefits accrue (μg/m³).
    upper_threshold       = Parameter()                      # Exposure threshold above which no health benefits accrue (μg/m³).
    β                     = Parameter()                      # Parameter linking a unit change in Pm2.5 exposure to a change in the risk of adult (>30) mortailty from all cases.
    SO2₀                  = Parameter(index=[regions])       # Initial value for baseline (pre-mitigation) SO₂ emissions.
    PM25₀                 = Parameter(index=[regions])       # Initial value for baseline (pre-mitigation) PM2.5 emissions.
    NOX₀                  = Parameter(index=[regions])       # Initial value for baseline (pre-mitigation) NOx emissions.
    exposure₀             = Parameter(index=[regions])       # Average exposure to PM2.5 in 2005 (μg/m³).
    source_receptor_SO2   = Parameter(index=[regions])       # Change in population-weighted PM2.5 exposure given a unit change in SO₂ emissions.
    source_receptor_PM25  = Parameter(index=[regions])       # Change in population-weighted PM2.5 exposure given a unit change in PM2.5 emissions.
    source_receptor_NOX   = Parameter(index=[regions])       # Change in population-weighted PM2.5 exposure given a unit change in NOx emissions.
    Θ                     = Parameter(index=[time, regions]) # Total years of life lost prematurely from all causes (millions/yr).
    base_deaths           = Parameter(index=[time, regions]) # Baseline deaths from all causes (millions/yr)
    base_emissions_SO2    = Parameter(index=[time, regions]) # Baseline (pre-mitigation) SO₂ emissions.
    base_emissions_PM25   = Parameter(index=[time, regions]) # Baseline (pre-mitigation) PM2.5 emissions.
    base_emissions_NOX    = Parameter(index=[time, regions]) # Baseline (pre-mitigation) NOx emissions.
    Δemissions_SO2        = Parameter(index=[time, regions]) # Change (co-reduction) in SO₂ emissions due to CO₂ mitigation.
    Δemissions_PM25       = Parameter(index=[time, regions]) # Change (co-reduction) in PM2.5 emissions due to CO₂ mitigation.
    Δemissions_NOX        = Parameter(index=[time, regions]) # Change (co-reduction) in NOx emissions due to CO₂ mitigation.

    ###############################################################################################################################
    # Model Variables.
    ###############################################################################################################################
    base_exposure         = Variable(index=[time, regions])  # Exposure to PM2.5 under baseline emissions (μg/m³).
    Δexposure             = Variable(index=[time, regions])  # Change in PM2.5 exposure due to CO₂ mitigation (μg/m³).
    lifeyears             = Variable(index=[time, regions])  # Life years gained from co-reduction of aerosols.
    avoided_deaths        = Variable(index=[time, regions])  # Avoided deaths from air pollution policy and co-reduction of aerosols (millions/yr)

    ###############################################################################################################################
    # Model Equations
    ###############################################################################################################################
    
    # Create run_timestep function to loop through each equation.
    function run_timestep(p, v, d, t)

        for r in d.regions

            # Calculate individual pollutant exposures as temporary variables, then save their sum as "base_exposure".
            exposure_SO2  = (p.SO2₀[r]  - p.base_emissions_SO2[t,r])  * p.source_receptor_SO2[r]
            exposure_PM25 = (p.PM25₀[r] - p.base_emissions_PM25[t,r]) * p.source_receptor_PM25[r]
            exposure_NOX  = (p.NOX₀[r]  - p.base_emissions_NOX[t,r])  * p.source_receptor_NOX[r]

            v.base_exposure[t,r] = max(p.exposure₀[r] - (exposure_SO2 + exposure_PM25 + exposure_NOX), 0.0)

            # Calculate change in pollutant exposures from CO₂ mitigation as temporary variables, then save their sum as "Δexposure".
            Δexposure_SO2  = p.Δemissions_SO2[t,r]  * p.source_receptor_SO2[r]
            Δexposure_PM25 = p.Δemissions_PM25[t,r] * p.source_receptor_PM25[r]
            Δexposure_NOX  = p.Δemissions_NOX[t,r]  * p.source_receptor_NOX[r]

            v.Δexposure[t,r] = Δexposure_SO2 + Δexposure_PM25 + Δexposure_NOX

            # Calculate life years, accounting for upper and lower PM2.5 exposure thresholds where health impacts accrue.

            # Find minimum between (i) base_exposure minus lower threshold and (ii) the change in exposure from CO₂ mitigation.
            A = min(v.base_exposure[t,r] - p.τ, v.Δexposure[t,r])

            # Find minimum between (i) upper and lower thresholds and (ii) A.
            B = min(p.upper_threshold - p.τ, A)

            # Find minimum between (i) exposure change - base exposure + upper threhold and (ii) B.
            C = min(v.Δexposure[t,r] - v.base_exposure[t,r] + p.upper_threshold, B)

            # Use C to calculate life years gained from co-reduction of aerosols.
            v.lifeyears[t,r] = p.β * p.Θ[t,r] * max(C, 0.0)

            # Use C to calculate avoided deaths from co-reduction of aerosols.
            v.avoided_deaths[t,r] = p.β * p.base_deaths[t,r] * max(C, 0.0)
        end
    end
end
