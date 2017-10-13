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
    mpcSol.cost = zeros(5)           # since in the Learning MPC we have 5 expressions composing the cost function, to be consistent with the dimension we initialize this with 5 elements...
    mpcSol.cost = [getvalue(m.costPath);getvalue(m.derivCost);getvalue(m.controlCost);0;0] # ... and set to 0 the last two elements
    
    nothing 
end   



function solveLearning_MPC(m::initLearningModel,mpcSol::classes.MpcSol,mpcParams::classes.MpcParams,trackCoeff::classes.TrackCoeff,modelParams::classes.ModelParams,zCurr::Array{Float64},uCurr::Array{Float64},selectedStates::classes.SelectedStates)

    println("FLAG LMPC")
 
    coeffCurvature  = trackCoeff.coeffCurvature::Array{Float64,1}
    selStates       = selectedStates.selStates::Array{Float64,2}
    statesCost      = selectedStates.statesCost::Array{Float64,1}
    sol_u::Array{Float64,2}    # where all the optimal control actions will be saved
    sol_z::Array{Float64,2}    # where all the optimal states will be saved

    # Set values needed by the MPC

    setvalue(m.z0,zCurr)
    println("value of Zcurr inside solveMPC= ",zCurr)
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
    mpcSol.cost  = zeros(5)           # since in the Learning MPC we have 5 expressions composing the cost function, to be consistent with the dimension we initialize this with 5 elements...
    mpcSol.cost  = [getvalue(m.derivCost);getvalue(m.controlCost);getvalue(m.laneCost);getvalue(m.slackCost);getvalue(m.terminalCost)] # ... and set to 0 the last two elements
    
    nothing 
end   