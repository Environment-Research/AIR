
@defcomp air_consumption begin

    ###############################################################################################################################
    # Model Indices
    ###############################################################################################################################
    regions   = Index()                                 # 12 regions in RICE2010.

    ###############################################################################################################################
    # Model Parameters
    ###############################################################################################################################
    use_VSL::Bool    = Parameter()                      # True/False indicator for whether or not to use a VSL with avoided deaths.
    Hyears           = Parameter()                      # Number of years of consumption to use in VOLY calculations.
    VOLY_elasticity  = Parameter()                      # Consumption elasticity of VOLY.
    pop              = Parameter(index=[time, regions]) # Regional population levels (millions of people).
    lifeyears        = Parameter(index=[time, regions]) # The number of life-years gained from CO2 mitigation (millions/yr).
    CPC              = Parameter(index=[time, regions]) # Per capita consumption ($1000s).
    avoided_deaths   = Parameter(index=[time, regions]) # Avoided deaths from co-reduction of aerosols (millions/yr).
    Y                = Parameter(index=[time, regions]) # Gross output minus climate damages and CO₂ mitigation costs).

    ###############################################################################################################################
    # Model Variables.
    ###############################################################################################################################
    VSL              = Variable(index=[time, regions])  # VSL ($1000s).
    health_benefit   = Variable(index=[time, regions])  # Total health benefit (billions $/yr).
    VOLY             = Variable(index=[time, regions])  # Value of selected life year type used to calculate monetized health benefits ($1000s/yr).
    CPC_with_health  = Variable(index=[time, regions])  # Per capita consumption levels with monetized health beneifits added ($1000s).

    ###############################################################################################################################
    # Model Equations
    ###############################################################################################################################
    
    # Create run_timestep function to loop through each equation.
    function run_timestep(p, v, d, t)

        for r in d.regions

            # Calculate the value of a life year based on regional consumption values.
            v.VOLY[t,r] = p.Hyears * p.CPC[t,r]^p.VOLY_elasticity

            if p.use_VSL
                # Calculate the monetized health benefit as VSL times avoided deaths.
                # Use VSL specification from Robinson, Hammit, & O'Keeffe (2017).
                v.VSL[t,r] = 9000.0 * ((p.Y[t,r]/p.pop[t,r]*1000) / 54.0)^1.4
                v.health_benefit[t,r] = p.avoided_deaths[t,r] * v.VSL[t,r] 
            else
                # Calculate the monetized health benefit as lifeyears gained times VOLY.
                v.health_benefit[t,r] = v.VOLY[t,r] * p.lifeyears[t,r]
            end
        
            # Calculate per capita consumption levels after adding monetized health benefits.
            v.CPC_with_health[t,r] = p.CPC[t,r] + (v.health_benefit[t,r] / p.pop[t,r])
        end
    end
end
