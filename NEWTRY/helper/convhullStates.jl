# This function selects, from the previous laps, the states that will be used to build the convex hull in the LMPC formulation. 
# Moreover, it returns their costs, as evaluated in the previous laps. 

# z[1] = s 
# z[2] = ey
# z[3] = epsi
# z[4] = v


function convhullStates(oldTraj::classes.OldTrajectory, posInfo::classes.PosInfo, mpcParams::classes.MpcParams,lapStatus::classes.LapStatus,selectedStates::classes.SelectedStates, modelParams::classes.ModelParams,simVariables::classes.SimulationVariables)

    # Read Inputs
    s               = posInfo.s          # current s 
    s_target        = posInfo.s_target   # s value of the target

    # Parameters
    Np              = selectedStates.Np  # number of states to choose per each last lap
    Nl              = selectedStates.Nl  # number of previous laps to take in consideration for the convex hull
    N               = mpcParams.N        # prediction horizon    
    dt              = modelParams.dt     # time step
    n_pf            = simVariables.n_pf
    
 
    selected_laps = zeros(Int64,Nl)
    for i = 1:Nl
        selected_laps[i] = lapStatus.currentLap-i    # use previous lap
    end
    
    #if lapStatus.currentLap >= 5
    #    selected_laps[2] = indmin(oldTraj.oldCost[2:lapStatus.currentLap-2])+1      # selects the best from all previous laps
    #end
 

    # Select the old data and the obstacle data
    oldS            = oldTraj.oldTraj[:,1,selected_laps]::Array{Float64,3}
    oldeY           = oldTraj.oldTraj[:,2,selected_laps]::Array{Float64,3}
    oldePsi         = oldTraj.oldTraj[:,3,selected_laps]::Array{Float64,3}
    oldV            = oldTraj.oldTraj[:,4,selected_laps]::Array{Float64,3}

    

    N_points        = size(oldTraj.oldTraj,1)     # first dimension = length = buffersize

    local s_total::Float64        # initialize
    local DistS::Array{Float64}   # initialize
    local idx_s::Array{Int64}     # initialize


    # Compute the total s (current position along track)
    s_total = s % s_target


    # Compute the index
    DistS = ( s_total - oldS ).^2   # compute all the distances between our current s and the s in the previous laps  

    idx_s = findmin(DistS,1)[2]     # find the indices of the nearest points in the old laps

    off = 5
    idx_s = idx_s + off

    
    
    for j = 0:(Nl-1)

        selectedStates.selStates[i=(j*Np)+1:(j+1)*Np,i=1:4] = oldTraj.oldTraj[i=idx_s[j+1]-(j*N_points):idx_s[j+1]+Np-(j*N_points)-1,i=1:4,selected_laps[j+1]]  # select the states from lap j...
        
        selectedStates.statesCost[i=(j*Np)+1:(j+1)*Np] = oldTraj.cost2target[i=idx_s[j+1]-(j*N_points):idx_s[j+1]-(j*N_points)+Np-1,selected_laps[j+1]]  # and their cost

        
        
        
    end



    nothing
end
