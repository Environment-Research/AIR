
@defcomp air_coreduction begin

    ###############################################################################################################################
    # Model Indices
    ###############################################################################################################################
    regions  = Index()                                          # 12 regions in RICE2010.
    aerosols = Index()                                          # 5 aerosols (SO₂, PM2.5, NOx, BC, OC).

    ###############################################################################################################################
    # Model Parameters
    ###############################################################################################################################
    MIU                  = Parameter(index=[time, regions])     # CO₂ abatement (share of total emissions reduced).
    base_emissions_SO2   = Parameter(index=[time, regions])     # Baseline (pre-mitigation) SO₂ emissions (Gg).
    base_emissions_PM25  = Parameter(index=[time, regions])     # Baseline (pre-mitigation) PM2.5 emissions (Gg).
    base_emissions_NOX   = Parameter(index=[time, regions])     # Baseline (pre-mitigation) NOx emissions (Gg).
    base_emissions_BC    = Parameter(index=[time, regions])     # Baseline (pre-mitigation) BC emissions (Gg).
    base_emissions_OC    = Parameter(index=[time, regions])     # Baseline (pre-mitigation) OC emissions (Gg).
    κ                    = Parameter(index=[regions, aerosols]) # Parameter describing effectiveness of CO₂ mitigation in co-reducing aerosol emissions (Ordered as: SO₂, PM2.5, NOx, BC, OC).

    ###############################################################################################################################
    # Model Variables.
    ###############################################################################################################################
    Δemissions_SO2       = Variable(index = [time, regions])    # Change (co-reduction) in SO₂ emissions due to CO₂ mitigation (Gg).
    Δemissions_PM25      = Variable(index = [time, regions])    # Change (co-reduction) in PM2.5 emissions due to CO₂ mitigation (Gg).
    Δemissions_NOX       = Variable(index = [time, regions])    # Change (co-reduction) in NOx emissions due to CO₂ mitigation (Gg).
    Δemissions_BC        = Variable(index = [time, regions])    # Change (co-reduction) in BC emissions due to CO₂ mitigation (Gg).
    Δemissions_OC        = Variable(index = [time, regions])    # Change (co-reduction) in OC emissions due to CO₂ mitigation (Gg).
    post_mitigation_SO2  = Variable(index = [time, regions])    # Post CO₂ mitigation emission levels of SO₂ (Gg).
    post_mitigation_PM25 = Variable(index = [time, regions])    # Post CO₂ mitigation emission levels of PM2.5 (Gg).
    post_mitigation_NOX  = Variable(index = [time, regions])    # Post CO₂ mitigation emission levels of NOx (Gg).
    post_mitigation_BC   = Variable(index = [time, regions])    # Post CO₂ mitigation emission levels of BC (Gg).
    post_mitigation_OC   = Variable(index = [time, regions])    # Post CO₂ mitigation emission levels of OC (Gg).

    ###############################################################################################################################
    # Model Equations
    ###############################################################################################################################
    
    # Create run_timestep function to loop through each equation.
    function run_timestep(p, v, d, t)

        for r in d.regions

            # Calculate change in aerosol emissions (co-reduction) due to CO₂ mitigation.
            v.Δemissions_SO2[t,r]  = p.κ[r,1] * p.base_emissions_SO2[t,r]  * p.MIU[t,r]
            v.Δemissions_PM25[t,r] = p.κ[r,2] * p.base_emissions_PM25[t,r] * p.MIU[t,r]
            v.Δemissions_NOX[t,r]  = p.κ[r,3] * p.base_emissions_NOX[t,r]  * p.MIU[t,r]
            v.Δemissions_BC[t,r]   = p.κ[r,4] * p.base_emissions_BC[t,r]   * p.MIU[t,r]
            v.Δemissions_OC[t,r]   = p.κ[r,5] * p.base_emissions_OC[t,r]   * p.MIU[t,r]

            # Calculate post CO₂ mitigation emission levels for each aerosol type.
            v.post_mitigation_SO2[t,r]  = p.base_emissions_SO2[t,r]  - v.Δemissions_SO2[t,r]
            v.post_mitigation_PM25[t,r] = p.base_emissions_PM25[t,r] - v.Δemissions_PM25[t,r]
            v.post_mitigation_NOX[t,r]  = p.base_emissions_NOX[t,r]  - v.Δemissions_NOX[t,r]
            v.post_mitigation_BC[t,r]   = p.base_emissions_BC[t,r]   - v.Δemissions_BC[t,r]
            v.post_mitigation_OC[t,r]   = p.base_emissions_OC[t,r]   - v.Δemissions_OC[t,r]
        end
    end
end
