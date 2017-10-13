# This file saves the data from the current lap in oldTrajectory



function saveOldTraj(oldTraj,cost_log::Array{Float64},zCurr::Array{Float64}, zCurr_x::Array{Float64},uCurr::Array{Float64},lapStatus::classes.LapStatus,simVariables::classes.SimulationVariables,z_pred_log::Array{Float64},u_pred_log::Array{Float64},modelParams::classes.ModelParams,curvature_curr::Array{Float64},alpha_log::Array{Float64},posInfo::classes.PosInfo)

    i               = lapStatus.currentIt           # current iteration number, just to make notation shorter
    costLap         = i                             # define as cost of the lap the number of iterations it took to arrive at the finish line
    currentLap      = lapStatus.currentLap   
    dt              = modelParams.dt    
    buffersize      = simVariables.buffersize 
    zCurr_export    = zeros(buffersize,4)
    uCurr_export    = zeros(buffersize,2)
    postbuff        = simVariables.postbuff
    s_target        = posInfo.s_target
    
    
    zCurr_export    = cat(1,zCurr[1:i,:], [zCurr[i,1]+collect(1:buffersize-i)*dt*zCurr[i,4] ones(buffersize-i,1)*zCurr[i,2:4]]) # extrapolate values for after the finish line so that the old trjectory has feasible values for the interpolation in the next round 
    uCurr_export    = cat(1,uCurr[1:i,:], zeros(buffersize-i,2)) 
    # if currentLap == 1
    #     zCurr_export    = cat(1,zCurr[1:i,:], [zCurr[i,1]+collect(1:buffersize-i)*dt*zCurr[i,4] ones(buffersize-i,1)*zCurr[i,2:4]]) # extrapolate values for after the finish line so that the old trjectory has feasible values for the interpolation in the next round 
    # else
    #     zCurr_export    = cat(1,zCurr[1:i,:], [zCurr[i,1]+collect(1:buffersize-i)*dt*zCurr[i,4] oldTraj.oldTraj[i+1:buffersize,2:4,currentLap-1]]) # extrapolate values for after the finish line so that the old trjectory has feasible values for the interpolation in the next round 
    # end
    # uCurr_export    = cat(1,uCurr[1:i,:], zeros(buffersize-i,2))


    cost2target     = zeros(buffersize) # array containing the cost to arrive from each point of the old trajectory to the target
    
    #save the terminal cost
    for j = 1:buffersize
        cost2target[j] = mpcParams.Q_cost*(costLap-j+1)  # why do i need Q_cost?
    end

        # if currentLap >1
        #     oldTraj.oldTraj[oldTraj.costLap[currentLap-1]+1:oldTraj.costLap[currentLap-1]+postbuff+1,2:4,currentLap-1] = zCurr_export[1:postbuff+1,2:4]
        #     oldTraj.oldTraj[oldTraj.costLap[currentLap-1]+1:oldTraj.costLap[currentLap-1]+postbuff+1,1,currentLap-1] = zCurr_export[1:postbuff+1,1] + s_target
        # end

    oldTraj.oldTraj[:,:,currentLap]                    = zCurr_export
    oldTraj.oldInput[:,:,currentLap]                   = uCurr_export
    oldTraj.oldTrajXY[1:size(zCurr_x)[1],:,currentLap] = zCurr_x
    oldTraj.oldNIter[currentLap]                       = costLap
    oldTraj.costs[:,:,currentLap]                      = cost_log[:,:,currentLap]
    oldTraj.z_pred_sol[:,:,:,currentLap]               = z_pred_log[:,:,:,currentLap]
    oldTraj.u_pred_sol[:,:,:,currentLap]               = u_pred_log[:,:,:,currentLap]
    oldTraj.cost2target[:,currentLap]                  = cost2target
    oldTraj.curvature[:,currentLap]                    = curvature_curr
    oldTraj.oldAlpha[:,:,currentLap]                   = alpha_log[:,:,currentLap]
    oldTraj.costLap[currentLap]                        = convert(Int64,costLap)            
       
    
end