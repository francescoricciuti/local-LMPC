# Variable definitions
# m.z_Ol[i,j] = z_OpenLoop, open loop prediction of the state, i = state, j = step

# States:
# i = 1 -> s
# i = 2 -> ey
# i = 3 -> epsi
# i = 4 -> v

function solvePathFollowMpc!(m::initPathFollowingModel,mpcSol::classes.MpcSol,mpcCoeff::classes.MpcCoeff,mpcParams::classes.MpcParams,trackCoeff::classes.TrackCoeff,
    lapStatus::classes.LapStatus,posInfo::classes.PosInfo,modelParams::classes.ModelParams,zCurr::Array{Float64},uCurr::Array{Float64}, iter::Int64)
    # Load Parameters
    coeffCurvature  = trackCoeff.coeffCurvature::Array{Float64,1}
    

    local sol_u::Array{Float64,2} 
    local sol_z::Array{Float64,2} 
    # Update current initial condition
    setvalue(m.z0,zCurr')

    # Update model values
    setvalue(m.coeff,coeffCurvature)
    setvalue(m.uCurr,uCurr)


    # Solve Problem and return solution
    sol_status  = solve(m.mdl)
    sol_u       = getvalue(m.u_Ol)
    sol_z       = getvalue(m.z_Ol)
    #  not sure i need this  mpcSol.eps[3,iter] = getvalue(m.eps)
    # c_print = getvalue(m.c)

    #if iter%20 ==0
    # println("curvature: $c_print")
    #println("Predicting until s = $(sol_z[end,1])") 
    #end
    
    # save data to solution class
    mpcSol.a_x = sol_u[1,1]
    mpcSol.d_f = sol_u[1,2]
    mpcSol.u   = sol_u
    mpcSol.z   = sol_z
    mpcSol.solverStatus = sol_status
    mpcSol.cost = zeros(5)
    mpcSol.cost = [getvalue(m.costPath);getvalue(m.derivCost);getvalue(m.controlCost);0;0]
    nothing #nothing to return
end



function solveLearningMpcProblem!(m::initLearningModel,mpcSol::classes.MpcSol,mpcCoeff::classes.MpcCoeff,mpcParams::classes.MpcParams,trackCoeff::classes.TrackCoeff,
    lapStatus::classes.LapStatus,posInfo::classes.PosInfo,modelParams::classes.ModelParams,zCurr::Array{Float64},uCurr::Array{Float64},iter::Int64,selectedStates::classes.SelectedStates)
    # Load Parameters
    coeffCurvature  = trackCoeff.coeffCurvature::Array{Float64,1}

    local sol_u::Array{Float64,2} 
    local sol_z::Array{Float64,2} 
    # Update current initial condition
    setvalue(m.z0,zCurr')

    # Update model values
    setvalue(m.coeff,coeffCurvature)
    setvalue(m.selStates,selectedStates.selStates)
    setvalue(m.statesCost,selectedStates.statesCost)
    setvalue(m.uCurr,uCurr)


    # Solve Problem and return solution
    sol_status  = solve(m.mdl)
    counter = 1
    while sol_status != :Optimal && counter <= 2
        sol_status  = solve(m.mdl)
        counter += 1
        println("Not solved optimally, trying again...")
        # println("state = ",zCurr)
    end


    sol_u       = getvalue(m.u_Ol)
    sol_z       = getvalue(m.z_Ol)
    #mpcSol.eps[:,iter] = getvalue(m.eps)
    # c_print = getvalue(m.c)

    #if iter%20 ==0
    # println("curvature: $c_print")
    #println("Predicting until s = $(sol_z[end,1])") 
    #end
    
    # safe data to solution class
    mpcSol.a_x = sol_u[1,1]
    mpcSol.d_f = sol_u[1,2]
    mpcSol.u   = sol_u
    mpcSol.z   = sol_z
    mpcSol.solverStatus = sol_status
    mpcSol.cost = zeros(5)
    mpcSol.cost = [getvalue(m.derivCost);getvalue(m.controlCost);getvalue(m.laneCost);getvalue(m.slackCost); getvalue(m.terminalCost)]
    if counter < 3
        return true
    else
        return false
    end
end

