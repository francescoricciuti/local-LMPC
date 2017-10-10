# This file Initializes the two MPC models (path following MPC and Learning MPC)



type initPathFollowingModel

    mdl::JuMP.Model

    z0::Array{JuMP.NonlinearParameter,1}
    coeff::Array{JuMP.NonlinearParameter,1}

    uCurr::Array{JuMP.NonlinearParameter,1}

    z_Ol::Array{JuMP.Variable,2}
    u_Ol::Array{JuMP.Variable,2}


    dsdt::Array{JuMP.NonlinearExpression,1}
    c::Array{JuMP.NonlinearExpression,1}

    derivCost::JuMP.NonlinearExpression
    costPath::JuMP.NonlinearExpression
    controlCost::JuMP.NonlinearExpression

    function initPathFollowingModel(mpcParams::classes.MpcParams,modelParams::classes.ModelParams,trackCoeff::classes.TrackCoeff)

        m = new()

        #### Initialize parameters
        s_target   = posInfo.s_target

        dt         = modelParams.dt                # time step
        L_a        = modelParams.l_A               # distance from CoM of the car to the front wheels
        L_b        = modelParams.l_B               # distance from CoM of the car to the rear wheels
        u_lb       = modelParams.u_lb              # lower bounds for the control inputs
        u_ub       = modelParams.u_ub              # upper bounds for the control inputs
        z_lb       = modelParams.z_lb              # lower bounds for the states
        z_ub       = modelParams.z_ub              # upper bounds for the states
        v_max      = modelParams.v_max 

        v_ref      = mpcParams.vPathFollowing              # reference velocity for the path following 
        N          = mpcParams.N                           # Prediction horizon
        QderivZ    = mpcParams.QderivZ::Array{Float64,1}   # weights for the derivative cost on the states
        QderivU    = mpcParams.QderivU::Array{Float64,1}   # weights for the derivative cost on the control inputs
        R          = mpcParams.R::Array{Float64,1}         # weights on the control inputs
        Q          = mpcParams.Q::Array{Float64,1}         # weights on the states for path following

        z_Init      = zeros(N+1,4)                   # define intial conditions
        z_Init[:,4] = 0.6*ones(N+1)

        
        ey_max      = trackCoeff.width/2           # bound for the state ey (distance from the center track). It is set as half of the width of the track for obvious reasons
        n_poly_curv = trackCoeff.nPolyCurvature    # polynomial degree for curvature approximation

        #### Defining reference trajectory for path following

        z_Ref       = cat(2,zeros(N+1,3),v_ref*ones(N+1,1))  # Reference trajectory: path following -> stay on line and keep constant velocity
        u_Ref       = zeros(N,2)

        #### Defining model, variables and constraints 

        mdl = Model(solver = IpoptSolver(print_level=0))

        @variable(mdl, z_Ol[i=1:(N+1),j=1:4])  # Define states (z=s,ey,epsi,v) 
        @variable(mdl, u_Ol[i=1:N,j=1:2])     # Define control inputs (u=a_x,d_f) 

        for i=1:2
            for j=1:N
                setlowerbound(u_Ol[j,i], u_lb[j,i])  # set upper...
                setupperbound(u_Ol[j,i], u_ub[j,i])  # ... and lower bound for control
            end
        end
        for i=1:4
            for j=1:N+1
                setlowerbound(z_Ol[j,i], z_lb[j,i])  # set upper...
                setupperbound(z_Ol[j,i], z_ub[j,i])  # ... and lower bound for states
            end
        end

        @NLparameter(mdl, uCurr[i=1:2] == 0)              # this will be filled with the input applied in the previous iteration, so that we can minimize the difference 
                                                          # between the input applied the last iteration and the first input applied this iteration
        @NLparameter(mdl, z0[i=1:4] == z_Init[1,i])       # initial conditions for the states
        @NLconstraint(mdl, [i=1:4], z_Ol[1,i] == z0[i])   # set initial conditions

       
        
        @NLparameter(mdl, coeff[i=1:n_poly_curv+1] == trackCoeff.coeffCurvature[i]) # define the coefficients needed to define track's curvature...
        @NLexpression(mdl, c[i = 1:N], sum{coeff[j]*z_Ol[i,1]^(n_poly_curv-j+1),j=1:n_poly_curv} + coeff[n_poly_curv+1]) # ... and define track's curvature
        
        # define expressions needed in the definition of the system's dynamics...
        @NLexpression(mdl, bta[i = 1:N],  atan( L_b / (L_a + L_b) * tan(u_Ol[i,2]) ) )         # beta -> car's slip angle 
        @NLexpression(mdl, dsdt[i = 1:N], z_Ol[i,4]*cos(z_Ol[i,3]+bta[i])/(1-z_Ol[i,2]*c[i]))  # time derivative of s

        # ... and define system's dynamics (seen as constraints by the optimization problem)

        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,1]  == z_Ol[i,1] + dt*dsdt[i]  )                                    # s
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,2]  == z_Ol[i,2] + dt*z_Ol[i,4]*sin(z_Ol[i,3]+bta[i])  )            # ey
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,3]  == z_Ol[i,3] + dt*(z_Ol[i,4]/L_a*sin(bta[i])-dsdt[i]*c[i])  )   # epsi
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,4]  == z_Ol[i,4] + dt*(u_Ol[i,1]))                                  # v
     
        #### Defining the cost function 

        # Derivative cost
        #----------------
        @NLexpression(mdl, derivCost, sum{QderivZ[j]*sum{(z_Ol[i+1,j]-z_Ol[i,j])^2,i=1:N},j=1:4} +
                                      sum{QderivU[j]*((uCurr[j]-u_Ol[1,j])^2+sum{(u_Ol[i+1,j]-u_Ol[i,j])^2,i=1:N-1}),j=1:2})

        # Control Input cost
        #-------------------
        @NLexpression(mdl, controlCost, 0.5*sum{R[j]*(sum{(u_Ol[i,j]-u_Ref[i,j])^2 ,i=1:N}) ,j=1:2})
 
        # Path following cost
        #--------------------
        @NLexpression(mdl, costPath, 0.5*sum{Q[i]*sum{(z_Ol[j,i]-z_Ref[j,i])^2,j=1:N+1},i=1:4})    

        # Overall Cost function (objective of the minimization)
        # -----------------------------------------------------
        @NLobjective(mdl, Min, costPath + derivCost + controlCost)


        #### Update model values

        m.mdl         = mdl
        m.coeff       = coeff       # curvature coefficients
        m.u_Ol        = u_Ol        # control inputs
        m.z_Ol        = z_Ol        # states
        m.uCurr       = uCurr       # initial conditions for the input
        m.z0          = z0          # initial conditions for the states
        m.derivCost   = derivCost   # derivative cost
        m.controlCost = controlCost # control cost
        m.costPath    = costPath    # path following cost

        sol_stat=solve(mdl)
        println("Finished solve 1: $sol_stat")
        sol_stat=solve(mdl)
        println("Finished solve 2: $sol_stat")

        return m
    end
end


# !!!!!!!!!!!!NEED TO DEFINE LEARNING MPC HERE!!!!!!!!!!!!!