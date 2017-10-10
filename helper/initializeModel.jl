type initPathFollowingModel
    mdl::JuMP.Model

    z0::Array{JuMP.NonlinearParameter,1}
    coeff::Array{JuMP.NonlinearParameter,1}
    uCurr::Array{JuMP.NonlinearParameter,1}

    z_Ol::Array{JuMP.Variable,2}
    u_Ol::Array{JuMP.Variable,2}


    dsdt::Array{JuMP.NonlinearExpression,1}
    c::Array{JuMP.NonlinearExpression,1}

    costPath::JuMP.NonlinearExpression
    derivCost::JuMP.NonlinearExpression
    controlCost::JuMP.NonlinearExpression





    function initPathFollowingModel(mpcParams::classes.MpcParams,modelParams::classes.ModelParams,trackCoeff::classes.TrackCoeff,mpcCoeff::classes.MpcCoeff,n_oldTraj::Int64)

        m = new()

        dt         = modelParams.dt      # time step
        L_a        = modelParams.l_A     # distance from CoM of the car to the front wheels
        L_b        = modelParams.l_B     # distance from CoM of the car to the rear wheels
        u_lb       = modelParams.u_lb    # lower bounds for the control inputs
        u_ub       = modelParams.u_ub    # upper bounds for the control inputs
        z_lb       = modelParams.z_lb    # lower bounds for the states
        z_ub       = modelParams.z_ub    # upper bounds for the states

 
        Q               = mpcParams.Q                   #Cost of states just for path following
	    QderivZ         = mpcParams.QderivZ::Array{Float64,1}
        QderivU         = mpcParams.QderivU::Array{Float64,1}        
        Q_velocity      = mpcParams.Q_velocity
        R               = mpcParams.R                   # weight for control, is always used but curently 0
        order           = mpcCoeff.order                # polynomial order of terminal constraints and cost approximation
        
        v_ref           = mpcParams.vPathFollowing      # reference velocity for the path following 

        N           = mpcParams.N                       # Prediction horizon
        ey_max      = trackCoeff.width/2                # bound for the state ey (distance from the center track). It is set as half of the width of the track for obvious reasons
        n_poly_curv = trackCoeff.nPolyCurvature         # polynomial degree for curvature approximation

        

        #### Create function-specific parameters

        s_target = 100 # the weight for s in pathfollowing is set to zero, so this value is not used
        local z_Ref::Array{Float64,2}
        z_Ref       = cat(2,s_target*ones(N+1,1),zeros(N+1,2),v_ref*ones(N+1,1))  # Reference trajectory: path following -> stay on line and keep constant velocity
        u_Ref       = zeros(N,2)
	    z_Init      = zeros(N+1,4)
        z_Init[:,4] = 0.6*ones(N+1)
	    v_max       = modelParams.v_max
        max_alpha   = modelParams.max_alpha


        #### Defining model, variables and constraints 

        mdl = Model(solver = IpoptSolver(print_level=0))#mu_strategy=adaptive,warm_start_init_point="yes"))#, max_cpu_time=0.08))#,linear_solver="ma57",max_iter=500, print_user_options="yes",max_cpu_time=2.0,))

        @variable(mdl, z_Ol[i=1:(N+1),j=1:4] ) # Define states (z=s,ey,epsi,v) with its upper and lower bounds
        @variable(mdl, u_Ol[i=1:N,j=1:2] ) # Define control inputs (u=a_x,d_f) with its upper and lower bounds

        for i=1:2
            for j=1:N
                setlowerbound(u_Ol[j,i], u_lb[j,i])
                setupperbound(u_Ol[j,i], u_ub[j,i])
            end
        end
        for i=1:4
            for j=1:N+1
                setlowerbound(z_Ol[j,i], z_lb[j,i])
                setupperbound(z_Ol[j,i], z_ub[j,i])
            end
        end

 
        @NLparameter(mdl, z0[i=1:4] == z_Init[1,i])       # initial conditions for the states
        @NLconstraint(mdl, [i=1:4], z_Ol[1,i] == z0[i])   # set the first state of the optimization as the initial conditions


        @NLparameter(mdl, coeff[i=1:n_poly_curv+1] == trackCoeff.coeffCurvature[i])
        @NLparameter(mdl, uCurr[i=1:2] == 0)

  
        
        @NLexpression(mdl, c[i = 1:N],    sum{coeff[j]*z_Ol[i,1]^(n_poly_curv-j+1),j=1:n_poly_curv} + coeff[n_poly_curv+1]) 

        @NLexpression(mdl, bta[i = 1:N],  atan( L_b / (L_a + L_b) * tan(u_Ol[i,2]) ) )         # beta -> car's slip angle 

        @NLexpression(mdl, dsdt[i = 1:N], z_Ol[i,4]*cos(z_Ol[i,3]+bta[i])/(1-z_Ol[i,2]*c[i]))  # time derivative of s

        # System dynamics (implemented as constraints)
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,1]  == z_Ol[i,1] + dt*dsdt[i]  )                                             # s
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,2]  == z_Ol[i,2] + dt*z_Ol[i,4]*sin(z_Ol[i,3]+bta[i])  )                     # ey
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,3]  == z_Ol[i,3] + dt*(z_Ol[i,4]/L_b*sin(bta[i])-dsdt[i]*c[i])  )            # epsi
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,4]  == z_Ol[i,4] + dt*(u_Ol[i,1]))#- 0.23*abs(z_Ol[i,4]) * z_Ol[i,4]))#0.63  # v
     


        # Cost function 

        # Derivative cost
        # ---------------------------------
        @NLexpression(mdl, derivCost, sum{QderivZ[j] * sum{(z_Ol[i,j] - z_Ol[i + 1,j]) ^ 2 , i = 1:N} ,j = 1:4} +
                                      sum{QderivU[j]*((uCurr[j]-u_Ol[1,j])^2+sum{(u_Ol[i-1,j]-u_Ol[i,j])^2 ,i=2:N}) ,j=1:2} )
       
        # Control Input cost
        # ---------------------------------
        @NLexpression(mdl, controlCost, 0.5*sum{R[j]*sum{(u_Ol[i,j]-u_Ref[i,j])^2 ,i=1:N} ,j=1:2})

        # ---------------------------------
        # Path Following cost
        # if we're in the first lap, just do path following
        @NLexpression(mdl, costPath, 0.5*sum{Q[i]*sum{(z_Ol[j,i]-z_Ref[j,i])^2,j=1:N+1} ,i=1:4})    # Follow trajectory

        # Overall Cost function (objective of the minimization)
        # ---------------------------------
        @NLobjective(mdl, Min, costPath + derivCost + controlCost)
        m.mdl = mdl

        # Update model values
        m.coeff = coeff # curvature coefficients
        m.uCurr = uCurr #last applied input
        m.u_Ol = u_Ol   #control inputs
        m.z_Ol = z_Ol   #states
        m.z0  = z0      #initial conditions
        m.derivCost= derivCost
        m.controlCost = controlCost 
        m.costPath = costPath
        

        return m
    end
end



type initLearningModel
    mdl::JuMP.Model

    #ssInfOn::Array{JuMP.NonlinearParameter,1}
    z0::Array{JuMP.NonlinearParameter,1}
    coeff::Array{JuMP.NonlinearParameter,1}
    uCurr::Array{JuMP.NonlinearParameter,1}
    #coeffTermConst::Array{JuMP.NonlinearParameter,3}
    #coeffTermCost::Array{JuMP.NonlinearParameter,2}
    selStates::Array{JuMP.NonlinearParameter,2}
    statesCost::Array{JuMP.NonlinearParameter,1}

    #s_startC::JuMP.NonlinearParameter

    z_Ol::Array{JuMP.Variable,2}
    u_Ol::Array{JuMP.Variable,2}
    #lambda::Array{JuMP.Variable,1}
    eps_sl::Array{JuMP.Variable,1}
    eps_alpha::Array{JuMP.Variable,1}
    alpha::Array{JuMP.Variable,1}

    dsdt::Array{JuMP.NonlinearExpression,1}
    c::Array{JuMP.NonlinearExpression,1}

    derivCost::JuMP.NonlinearExpression
    controlCost::JuMP.NonlinearExpression
    terminalCost::JuMP.NonlinearExpression

    #soft constraints

    laneCost::JuMP.NonlinearExpression
    slackCost::JuMP.NonlinearExpression


    function initLearningModel(mpcParams::classes.MpcParams,modelParams::classes.ModelParams,trackCoeff::classes.TrackCoeff,mpcCoeff::classes.MpcCoeff, n_oldTraj,selectedStates::classes.SelectedStates)

        m = new()

        dt        = modelParams.dt
        L_a       = modelParams.l_A
        L_b       = modelParams.l_B
        u_lb      = modelParams.u_lb
        u_ub      = modelParams.u_ub
        z_lb      = modelParams.z_lb
        z_ub      = modelParams.z_ub

        Q               = mpcParams.Q                #Cost of states just for path following
        Q_term          = mpcParams.Q_term
        Q_cost          = mpcParams.Q_cost
        Q_lane          = mpcParams.Q_lane
        Q_velocity      = mpcParams.Q_velocity
        R               = mpcParams.R                # cost for control is always used but curently 0
        
        selStates       = selectedStates.selStates::Array{Float64,2}
        statesCost      = selectedStates.statesCost::Array{Float64}
        Np              = selectedStates.Np::Int64
        
        order           = mpcCoeff.order       # polynomial order of terminal constraints and cost approximation
        QderivZ         = mpcParams.QderivZ::Array{Float64,1}
        QderivU         = mpcParams.QderivU::Array{Float64,1}
        v_ref           = mpcParams.vPathFollowing

        N           = mpcParams.N
        ey_max      = trackCoeff.width/2
        n_poly_curv = trackCoeff.nPolyCurvature         # polynomial degree of curvature approximation

        #### Create function-specific parameters

        u_Ref           = zeros(N,2)

        z_Init          = zeros(N+1,4)
        z_Init[:,4]     = v_ref*ones(N+1)

        v_max = modelParams.v_max
        max_alpha =modelParams.max_alpha


        #### Defining model, variables and constraints

        mdl = Model(solver = IpoptSolver(print_level=0))#mu_strategy=adaptive,warm_start_init_point="yes"))#, max_cpu_time=0.08))#,linear_solver="ma57",max_iter=500, print_user_options="yes",max_cpu_time=2.0,))

        @variable(mdl,z_Ol[i=1:(N+1),j=1:4] ) # Define states (z=s,ey,epsi,v) with its upper and lower bounds
        @variable(mdl,u_Ol[i=1:N,j=1:2] ) # Define control inputs (u=a_x,d_f) with its upper and lower bounds
        @variable(mdl, alpha[1:2*Np] <= 1)                                                 # coefficients of the convex hull
        @variable(mdl, eps_sl[1:N+1] >= 0)                                                         # eps for soft lane constraints
        @variable(mdl, eps_alpha[1:4] >=0)                                                      # eps for soft constraint on alpha

        for i=1:2
            for j=1:N
                setlowerbound(u_Ol[j,i], u_lb[j,i])
                setupperbound(u_Ol[j,i], u_ub[j,i])
            end
        end
        for i=1:4
            for j=1:N+1
                setlowerbound(z_Ol[j,i], z_lb[j,i])
                setupperbound(z_Ol[j,i], z_ub[j,i])
            end
        end

        @NLparameter(mdl, selStates[1:2*Np,1:4] == 0)     # states from the previous trajectories selected in "convhullStates"
        @NLparameter(mdl, statesCost[1:2*Np] == 0)        # costs of the states selected in "convhullStates"
 
        @NLparameter(mdl, z0[i=1:4] == z_Init[1,i])       # initial conditions for the states

        @NLconstraint(mdl, [i=1:4], z_Ol[1,i] == z0[i])   # set the first state of the optimization as the initial conditions

        #@NLparameter(mdl,ssInfOn[1:n_oldTraj]== 1)


        @NLconstraint(mdl, sum{alpha[i],i=1:2*Np} == 1)    # constraint on the coefficients of the convex hull

        for i = 1:4
            #@NLconstraint(mdl,z_Ol[N+1,i] == sum{alpha[j]*selStates[j,i],j=1:numb_laps*Np})
            @NLconstraint(mdl,z_Ol[N+1,i] >= sum{alpha[j]*selStates[j,i],j=1:2*Np}-eps_alpha[i])  # terminal constraint are implemented as soft constraints
            @NLconstraint(mdl,z_Ol[N+1,i] <= sum{alpha[j]*selStates[j,i],j=1:2*Np}+eps_alpha[i])  #  terminal constraint are implemented as soft constraints
        end     


        @NLconstraint(mdl, [i=1:N+1], z_Ol[i,2] <= ey_max + eps_sl[i])   # lane constraint are implemented as soft constraints 
        @NLconstraint(mdl, [i=1:N+1], z_Ol[i,2] >= -ey_max - eps_sl[i])  # lane constraint are implemented as soft constraints 


        @NLparameter(mdl, coeff[i=1:n_poly_curv+1] == trackCoeff.coeffCurvature[i])
        @NLparameter(mdl, uCurr[i=1:2] == 0)
       
        @NLexpression(mdl, c[i = 1:N],    coeff[1]*z_Ol[i,1]^4+coeff[2]*z_Ol[i,1]^3+coeff[3]*z_Ol[i,1]^2+coeff[4]*z_Ol[i,1]+coeff[5])
        #@NLexpression(mdl, c[i = 1:N],    sum{coeff[j]*z_Ol[i,1]^(n_poly_curv-j+1),j=1:n_poly_curv} + coeff[n_poly_curv+1]) 

        @NLexpression(mdl, bta[i = 1:N],  atan( L_b / (L_a + L_b) * tan(u_Ol[i,2]) ) )          # beta -> car's slip angle
        @NLexpression(mdl, dsdt[i = 1:N], z_Ol[i,4]*cos(z_Ol[i,3]+bta[i])/(1-z_Ol[i,2]*c[i]))   # time derivative of s


        # System dynamics (implemented as constraints)
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,1]  == z_Ol[i,1] + dt*dsdt[i]  )                                             # s
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,2]  == z_Ol[i,2] + dt*z_Ol[i,4]*sin(z_Ol[i,3]+bta[i])  )                     # ey
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,3]  == z_Ol[i,3] + dt*(z_Ol[i,4]/L_b*sin(bta[i])-dsdt[i]*c[i])  )            # epsi
        @NLconstraint(mdl, [i=1:N], z_Ol[i+1,4]  == z_Ol[i,4] + dt*(u_Ol[i,1]))#- 0.23*abs(z_Ol[i,4]) * z_Ol[i,4]))#0.63  # v
     
        # @NLconstraint(mdl,[i =1:(N+1)],( (z_Ol[i,1]-sCoord_obst[i,1])/rs )^2 + ( (z_Ol[i,2]-sCoord_obst[i,2])/ry) ^2 - 1>=0)

        # Cost function 

        # Derivative cost
        # ---------------------------------
        @NLexpression(mdl, derivCost, sum{QderivZ[j]*(sum{(z_Ol[i,j]-z_Ol[i+1,j])^2,i=1:N}),j=1:4} )#+
                                          #QderivU[1]*( ((uCurr[1]-u_Ol[1,1])^2)+sum{(u_Ol[i-1,1]-u_Ol[i,1])^2,i=2:(N)})+
                                          #QderivU[2]*( ((uCurr[2]-u_Ol[1,2])^2)+sum{(u_Ol[i-1,2]-u_Ol[i,2])^2,i=2:(N)}))
       # @NLexpression(mdl, derivCost, QderivU[1]*((uCurr[1]-u_Ol[1,1])^2+sum{(u_Ol[j,1]-u_Ol[j+1,1])^2, j=1:Ns}))

        # Lane cost
        # ---------------------------------
        #@NLexpression(mdl, laneCost, Q_lane*sum{z_Ol[i,2]^2*((0.5+0.5*tanh(35*(z_Ol[i,2]-ey_max-0.09))) + (0.5-0.5*tanh(35*(z_Ol[i,2]+ey_max+0.09)))),i=1:N+1})
        @NLexpression(mdl, laneCost, Q_lane*sum{10.0*eps_sl[i]+100.0*eps_sl[i]^2 ,i=1:2})


        # Control Input cost
        # --------------------------------- 
        @NLexpression(mdl, controlCost, 0.5*sum{R[j]*sum{(u_Ol[i,j]-u_Ref[i,j])^2 ,i=1:N} ,j=1:2})

        # Slack cost (soft)
        # ---------------------------------
        @NLexpression(mdl, slackCost, sum{10*eps_alpha[i]+100*eps_alpha[i]^2,i=1:4})

        # Terminal Cost
        # ---------------------------------
        @NLexpression(mdl, terminalCost , sum{alpha[i]*statesCost[i], i=1:2*Np})


        # State cost
        # ---------------------------------
        #@NLexpression(mdl, costZ, Q_cost*sum(1 for i=1:N+1))

        
        @NLobjective(mdl, Min, derivCost + laneCost + controlCost + slackCost + terminalCost)


        m.mdl = mdl
        m.coeff = coeff # curvature coefficients
        m.selStates = selStates
        m.statesCost = statesCost
        m.uCurr = uCurr #last applied input
        m.u_Ol = u_Ol
        m.z_Ol = z_Ol
        #m.ssInfOn = ssInfOn
        m.eps_sl = eps_sl
        m-eps_alpha = eps_alpha
        m.z0  = z0

        m.derivCost = derivCost
        m.laneCost = laneCost
        m.controlCost = controlCost
        m.slackCost = slackCost
        m.terminalCost = terminalCost



        return m
    end
end
