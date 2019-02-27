
@defcomp welfare begin

    ###############################################################################################################################
    # Model Indices
    ###############################################################################################################################
    regions  = Index()                                 # 12 regions in RICE2010.

    ###############################################################################################################################
    # Model Parameters
    ###############################################################################################################################
    rho             = Parameter()                      # Pure rate of time preference.
    eta             = Parameter()                      # Elasticity of marginal utility of consumption.
    pop             = Parameter(index=[time, regions]) # Regional population levels (millions of people).
    CPC_with_health = Parameter(index=[time, regions]) # Per capita consumption levels with monetized health beneifits added ($1000s).

    ###############################################################################################################################
    # Model Variables.
    ###############################################################################################################################
    welfare         = Variable()                       # Total welfare over model time horizon.
    
    ###############################################################################################################################
    # Model Equations
    ###############################################################################################################################
    
    # Create run_timestep function to loop through each equation.
    function run_timestep(p, v, d, t)

        # Calculate welfare for each period.  Note, if eta = 1, the social welfare function becomes log(x).
        if is_first(t)
            if p.eta == 1.0
                v.welfare = sum(log(p.CPC_with_health[t,:]) .* p.pop[t,:]) / (1.0 + p.rho)^(10*(t.t - 1))
            else
                v.welfare = sum((p.CPC_with_health[t,:] .^ (1.0 - p.eta)) ./ (1.0 - p.eta) .* p.pop[t,:]) / (1.0 + p.rho)^(10*(t.t - 1))
            end
        else
            if p.eta == 1.0
                v.welfare = v.welfare + sum(log(p.CPC_with_health[t,:]) .* p.pop[t,:]) / (1.0 + p.rho)^(10*(t.t - 1))
            else
                v.welfare = v.welfare + sum((p.CPC_with_health[t,:] .^ (1.0 - p.eta)) ./ (1.0 - p.eta) .* p.pop[t,:]) / (1.0 + p.rho)^(10*(t.t - 1))
            end
        end
    end
end
