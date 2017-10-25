# This file solves the two MPC models defined in "initializeModels.jl"


function solvePF_MPC(m::initPathFollowingModel,mpcSol::classes.MpcSol,mpcParams::classes.MpcParams,trackCoeff::classes.TrackCoeff,modelParams::classes.ModelParams,zCurr::Array{Float64},uCurr::Array{Float64})

    # Define needed variables

    coeffCurvature  = trackCoeff.coeffCurvature::Array{Float64,1}
    sol_u::Array{Float64,2}    # where all the optimal control actions will be saved
    sol_z::Array{Float64,2}    # where all the optimal states will be saved

    # Set values needed by the MPC

    setvalue(m.z0,zCurr)
    setvalue(m.coeff,coeffCurvature)
    setvalue(m.uCurr,uCurr)


    # Solve Problem and return solution
    sol_status  = solve(m.mdl)
    sol_u       = getvalue(m.u_Ol)
    sol_z       = getvalue(m.z_Ol)
  
    
    # save data to solution class

    mpcSol.a_x = sol_u[1,1]         # we have a control horizon = 1, then we apply just the first input
    mpcSol.d_f = sol_u[1,2]         # same as above
    mpcSol.u   = sol_u              # here we save all the optimal inputs computed by the MPC
    mpcSol.z   = sol_z              # here we save all the optimal states computed by the MPC
    mpcSol.solverStatus = sol_status
    mpcSol.cost = zeros(7)           # since in the Learning MPC we have 7 expressions composing the cost function, to be consistent with the dimension we initialize this with 7 elements...
    mpcSol.cost = [getvalue(m.costPath);getvalue(m.derivCost);getvalue(m.controlCost);0;0;0;0] # ... and set to 0 the last four elements
    
    nothing 
end   



function solveLearning_MPC(m::initLearningModel,mpcSol::classes.MpcSol,mpcParams::classes.MpcParams,trackCoeff::classes.TrackCoeff,modelParams::classes.ModelParams,zCurr::Array{Float64},uCurr::Array{Float64},selectedStates::classes.SelectedStates)

    
 
    coeffCurvature  = trackCoeff.coeffCurvature::Array{Float64,1}
    selStates       = selectedStates.selStates::Array{Float64,2}
    statesCost      = selectedStates.statesCost::Array{Float64,1}
    sol_u::Array{Float64,2}    # where all the optimal control actions will be saved
    sol_z::Array{Float64,2}    # where all the optimal states will be saved

    # Set values needed by the MPC

    setvalue(m.z0,zCurr)
    setvalue(m.coeff,coeffCurvature)
    setvalue(m.uCurr,uCurr)
    setvalue(m.selStates,selStates)
    setvalue(m.statesCost,statesCost)

    # Solve Problem and return solution
    sol_status  = solve(m.mdl)
    sol_u       = getvalue(m.u_Ol)
    sol_z       = getvalue(m.z_Ol)
    sol_alpha   = getvalue(m.alpha)

    # save data to solution class

    mpcSol.a_x   = sol_u[1,1]         # we have a control horizon = 1, then we apply just the first input
    mpcSol.d_f   = sol_u[1,2]         # same as above
    mpcSol.u     = sol_u              # here we save all the optimal inputs computed by the MPC
    mpcSol.z     = sol_z              # here we save all the optimal states computed by the MPC
    mpcSol.alpha = sol_alpha
    mpcSol.solverStatus = sol_status
    mpcSol.cost  = zeros(7)           # since in the Obstacle Learning MPC we have 7 expressions composing the cost function, to be consistent with the dimension we initialize this with 7 elements...
    mpcSol.cost  = [getvalue(m.derivCost);getvalue(m.controlCost);getvalue(m.laneCost);getvalue(m.slackCost);getvalue(m.terminalCost);0;0] # ... and set to 0 the last two elements
    
    nothing 
end   



function solveObs_LMPC(m::initObsModel,mpcSol::classes.MpcSol,mpcParams::classes.MpcParams,trackCoeff::classes.TrackCoeff,modelParams::classes.ModelParams,zCurr::Array{Float64},uCurr::Array{Float64},selectedStates::classes.SelectedStates,obs_now::Array{Float64},obstacle::classes.Obstacle)

    
 
    coeffCurvature  = trackCoeff.coeffCurvature::Array{Float64,1}
    selStates       = selectedStates.selStates::Array{Float64,2}
    statesCost      = selectedStates.statesCost::Array{Float64,1}
    Q_obs           = mpcParams.Q_obs::Array{Float64,1}
    N               = mpcParams.N::Int64
    sol_u::Array{Float64,2}    # where all the optimal control actions will be saved
    sol_z::Array{Float64,2}    # where all the optimal states will be saved

    obs      = zeros(N+1,3)
    obs[1,:] = obs_now

    # Compute the position of the obstacle  in the whole prediction horizon

    for i = 1:N
        obs[i+1,1] = obs[i,1] + dt*i*obs[i,3]
        obs[i+1,2] = obs[i,2]
        obs[i+1,3] = obs[i,3]
    end

     #println("obs in solveMPC= ",obs)

    # Set values needed by the MPC

    setvalue(m.z0,zCurr)
    setvalue(m.coeff,coeffCurvature)
    setvalue(m.uCurr,uCurr)
    setvalue(m.selStates,selStates)
    setvalue(m.statesCost,statesCost)
    setvalue(m.Q_obs,Q_obs)
    setvalue(m.obs,obs)

    # Solve Problem and return solution
    sol_status  = solve(m.mdl)
    sol_u       = getvalue(m.u_Ol)
    sol_z       = getvalue(m.z_Ol)
    sol_alpha   = getvalue(m.alpha)

   # println("optimal alpha= ", sol_alpha)

    # save data to solution class

    mpcSol.a_x   = sol_u[1,1]         # we have a control horizon = 1, then we apply just the first input
    mpcSol.d_f   = sol_u[1,2]         # same as above
    mpcSol.u     = sol_u              # here we save all the optimal inputs computed by the MPC
    mpcSol.z     = sol_z              # here we save all the optimal states computed by the MPC
    mpcSol.alpha = sol_alpha
    mpcSol.solverStatus = sol_status
    mpcSol.cost  = zeros(7)           
    mpcSol.cost  = [getvalue(m.derivCost);getvalue(m.controlCost);getvalue(m.laneCost);0;getvalue(m.terminalCost);getvalue(m.obstacleCost1);getvalue(m.obstacleCost2)]
    println("obstacleCost1= ",getvalue(m.obstacleCost1))
    println("obstacleCost2= ",getvalue(m.obstacleCost2))
    nothing 
end   