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
        println("Finished solve 1 of PF_MPC: $sol_stat")
        sol_stat=solve(mdl)
        println("Finished solve 2 of PF_MPC: $sol_stat")

        return m
    end
end


# !!!!!!!!!!!!NEED TO DEFINE LEARNING MPC HERE!!!!!!!!!!!!!


type initLearningModel

    mdl::JuMP.Model

    z0::Array{JuMP.NonlinearParameter,1}
    coeff::Array{JuMP.NonlinearParameter,1}
    selStates::Array{JuMP.NonlinearParameter,2}
    statesCost::Array{JuMP.NonlinearParameter,1}


    eps_lane::Array{JuMP.Variable,1}
    eps_alpha::Array{JuMP.Variable,1}
    alpha::Array{JuMP.Variable,1}
    z_Ol::Array{JuMP.Variable,2}
    u_Ol::Array{JuMP.Variable,2}

    uCurr::Array{JuMP.NonlinearParameter,1}

    dsdt::Array{JuMP.NonlinearExpression,1}
    c::Array{JuMP.NonlinearExpression,1}

    derivCost::JuMP.NonlinearExpression
    controlCost::JuMP.NonlinearExpression
    slackCost::JuMP.NonlinearExpression
    laneCost::JuMP.NonlinearExpression
    terminalCost::JuMP.NonlinearExpression

    function initLearningModel(mpcParams::classes.MpcParams,modelParams::classes.ModelParams,trackCoeff::classes.TrackCoeff,selectedStates::classes.SelectedStates)

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

        ey_max      = trackCoeff.width/2           # bound for the state ey (distance from the center track). It is set as half of the width of the track for obvious reasons
        n_poly_curv = trackCoeff.nPolyCurvature    # polynomial degree for curvature approximation
        v_max       = modelParams.v_max 

        v_ref      = mpcParams.vPathFollowing              # reference velocity for the path following 
        N          = mpcParams.N                           # Prediction horizon
        QderivZ    = mpcParams.QderivZ::Array{Float64,1}   # weights for the derivative cost on the states
        QderivU    = mpcParams.QderivU::Array{Float64,1}   # weights for the derivative cost on the control inputs
        R          = mpcParams.R::Array{Float64,1}         # weights on the control inputs
        Q          = mpcParams.Q::Array{Float64,1}         # weights on the states for path following
        Q_lane     = mpcParams.Q_lane::Float64             # weight on the soft constraint on the lane
        Q_alpha    = mpcParams.Q_alpha::Float64            # weight on the soft constraint for the convex hull

      
        Np         = selectedStates.Np::Int64              # how many states to select


        #### Create function-specific parameters

        u_Ref           = zeros(N,2)

        z_Init          = zeros(N+1,4)
        z_Init[:,4]     = v_ref*ones(N+1)


        #### Defining model, variables and constraints

        mdl = Model(solver = IpoptSolver(print_level=0))#mu_strategy=adaptive,warm_start_init_point="yes"))#, max_cpu_time=0.08))#,linear_solver="ma57",max_iter=500, print_user_options="yes",max_cpu_time=2.0,))

        @variable(mdl, z_Ol[i=1:(N+1),j=1:4] ) # Define states (z=s,ey,epsi,v) with its upper and lower bounds
        @variable(mdl, u_Ol[i=1:N,j=1:2] )     # Define control inputs (u=a_x,d_f) with its upper and lower bounds
        @variable(mdl, alpha[1:2*Np] <= 1)     # coefficients of the convex hull
        @variable(mdl, eps_lane[1:N+1] >= 0)   # eps for soft lane constraints
        @variable(mdl, eps_alpha[1:4] >=0)     # eps for soft constraint on alpha


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


        @NLparameter(mdl, selStates[1:2*Np,1:4] == 0)                                 # states from the previous trajectories selected in "convhullStates"
        @NLparameter(mdl, statesCost[1:2*Np] == 0)                                    # costs of the states selected in "convhullStates"
        @NLparameter(mdl, z0[i=1:4] == z_Init[1,i])                                   # initial conditions for the states
        @NLconstraint(mdl, [i=1:4], z_Ol[1,i] == z0[i])                               # set initial conditions
        @NLparameter(mdl, coeff[i=1:n_poly_curv+1] == trackCoeff.coeffCurvature[i])   # coefficients for the curvature
        @NLparameter(mdl, uCurr[i=1:2] == 0) 
                                               # initial conditions for the control actions


        @NLconstraint(mdl, sum{alpha[i],i=1:2*Np} == 1)    # constraint on the coefficients of the convex hull
        for i = 1:4
            @NLconstraint(mdl,z_Ol[N+1,i] == sum{alpha[j]*selStates[j,i],j=1:2*Np})
            #@NLconstraint(mdl,z_Ol[N+1,i] >= sum{alpha[j]*selStates[j,i],j=1:2*Np}-eps_alpha[i])  # terminal constraint are implemented as soft constraints
            #@NLconstraint(mdl,z_Ol[N+1,i] <= sum{alpha[j]*selStates[j,i],j=1:2*Np}+eps_alpha[i])  #  terminal constraint are implemented as soft constraints
        end     
        @NLconstraint(mdl, [i=1:N+1], z_Ol[i,2] <= ey_max + eps_lane[i])   # lane constraint are implemented as soft constraints 
        @NLconstraint(mdl, [i=1:N+1], z_Ol[i,2] >= -ey_max - eps_lane[i])  # lane constraint are implemented as soft constraints 


        
        # define expressions needed in the definition of the system's dynamics...
        @NLexpression(mdl, c[i = 1:N], sum{coeff[j]*z_Ol[i,1]^(n_poly_curv-j+1),j=1:n_poly_curv} + coeff[n_poly_curv+1]) # track's curvature
        @NLexpression(mdl, bta[i = 1:N],  atan( L_b / (L_a + L_b) * tan(u_Ol[i,2]) ) )                                   # beta -> car's slip angle 
        @NLexpression(mdl, dsdt[i = 1:N], z_Ol[i,4]*cos(z_Ol[i,3]+bta[i])/(1-z_Ol[i,2]*c[i]))                            # time derivative of s

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

        # Slack cost (soft)
        # ---------------------------------
        @NLexpression(mdl, slackCost, Q_alpha*sum{10*eps_alpha[i]+100*eps_alpha[i]^2,i=1:4})

        # Lane cost (soft)
        # ---------------------------------
        @NLexpression(mdl, laneCost, Q_lane*sum{10.0*eps_lane[i]+100.0*eps_lane[i]^2 ,i=1:2})

        # Terminal Cost
        # ---------------------------------
        @NLexpression(mdl, terminalCost , sum{alpha[i]*statesCost[i], i=1:2*Np})

        # Overall Cost function (objective of the minimization)
        # -----------------------------------------------------

        @NLobjective(mdl, Min, derivCost + laneCost + controlCost + terminalCost )#+ slackCost)

        #### Update model values

        m.mdl         = mdl
        m.coeff       = coeff       # curvature coefficients
        m.u_Ol        = u_Ol        # control inputs
        m.z_Ol        = z_Ol        # states
        m.uCurr       = uCurr       # initial conditions for the input
        m.z0          = z0          # initial conditions for the states
        m.derivCost   = derivCost   # derivative cost
        m.controlCost = controlCost # control cost
        m.laneCost    = laneCost    # lane cost
        m.slackCost   = slackCost   # cost on alpha
        m.terminalCost= terminalCost# terminal cost
        m.selStates   = selStates   # selected states
        m.statesCost  = statesCost  # cost of the selected states
        m.alpha       = alpha       # parameters of the convex hull


        sol_stat=solve(mdl)
        println("Finished solve 1 of LMPC: $sol_stat")
        sol_stat=solve(mdl)
        println("Finished solve 2 of LMPC: $sol_stat")

        return m
    end
end
