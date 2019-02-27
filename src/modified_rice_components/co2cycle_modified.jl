
@defcomp co2cycle begin
   
    ###############################################################################################################################
    # Model Parameters
    ###############################################################################################################################
    mat0   = Parameter()             # Initial carbon concentration in atmosphere (GtC).
    mu0    = Parameter()             # Initial carbon concentration in upper strata (GtC).
    ml0    = Parameter()             # Initial carbon concentration in lower strata (GtC).
    b12    = Parameter()             # Carbon cycle transition matrix.
    b23    = Parameter()             # Carbon cycle transition matrix.
    b11    = Parameter()             # Carbon cycle transition matrix.
    b21    = Parameter()             # Carbon cycle transition matrix.
    b22    = Parameter()             # Carbon cycle transition matrix.
    b32    = Parameter()             # Carbon cycle transition matrix.
    b33    = Parameter()             # Carbon cycle transition matrix.
    E      = Parameter(index=[time]) # Total carbon emissions (GtC/year).

    ###############################################################################################################################
    # Model Variables
    ###############################################################################################################################
    MAT    = Variable(index=[time])  # Carbon concentration increase in atmosphere (GtC).
    MATSUM = Variable(index=[time])  # Sum of MAT[t] and MAT[t+1] to use in FORC[t] for radiativeforcing component.
    MU     = Variable(index=[time])  # Carbon concentration increase in shallow oceans (GtC).
    ML     = Variable(index=[time])  # Carbon concentration increase in lower oceans (GtC).

    ###############################################################################################################################
    # Model Equations
    ###############################################################################################################################
    
    # Create run_timestep function to loop through each equation.
    function run_timestep(p, v, d, t)
        
        #Define initial values and equations for for MAT, ML, MU, and MATSUM.
        if is_first(t)
            v.MAT[t] = p.mat0
            v.ML[t] = p.ml0
            v.MU[t] = p.mu0
            v.MATSUM[t] = 0
        else
            v.MAT[t] = v.MAT[t-1] * p.b11 + v.MU[t-1] * p.b21 + (p.E[t-1] * 10)
            v.ML[t] = v.ML[t-1] * p.b33 + v.MU[t-1] * p.b23
            v.MU[t] = v.MAT[t-1] * p.b12 + v.MU[t-1] * p.b22 + v.ML[t-1] * p.b32
            v.MATSUM[t] = v.MAT[t] + (v.MAT[t] * p.b11 + v.MU[t] * p.b21 +  (p.E[t] * 10))
        end
    end
end
