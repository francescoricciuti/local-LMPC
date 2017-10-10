# This function selects, from the previous laps, the states that will be used to build the convex hull in the LMPC formulation. 
# Moreover, it returns their costs, as evaluated in the previous laps. 

# structure of oldTrajectory: 1st dimension = step number (time equiv.), 2nd dimension = state number, 3rd dimennsion = lap numbe

# z[1] = s 
# z[2] = ey
# z[3] = epsi
# z[4] = v


function convhullStates(oldTraj::classes.OldTrajectory, mpcCoeff::classes.MpcCoeff, posInfo::classes.PosInfo, mpcParams::classes.MpcParams,lapStatus::classes.LapStatus,selectedStates::classes.SelectedStates)

    # Read Inputs
    s_start         = posInfo.s_start
    s               = posInfo.s
    s_target        = posInfo.s_target

    # Parameters
    Np              = selectedStates.Np      # number of states to choose per each last lap
    N               = mpcParams.N
 

    selected_laps = zeros(Int64,2)
    selected_laps[1] = lapStatus.currentLap-1                                   # use previous lap
    selected_laps[2] = lapStatus.currentLap-2                                   # and the one before
    
    #if lapStatus.currentLap >= 5
    #    selected_laps[2] = indmin(oldTraj.oldCost[2:lapStatus.currentLap-2])+1      # selects the best from all previous laps
    #end

    

    # Select the old data
    oldS            = oldTraj.oldTraj[:,1,selected_laps]::Array{Float64,3}
    oldeY           = oldTraj.oldTraj[:,2,selected_laps]::Array{Float64,3}
    oldePsi         = oldTraj.oldTraj[:,3,selected_laps]::Array{Float64,3}
    oldV            = oldTraj.oldTraj[:,4,selected_laps]::Array{Float64,3}


    N_points        = size(oldTraj.oldTraj,1)     # first dimension = length (=buffersize)

    local s_total::Float64        # initialize
    local DistS::Array{Float64}   # initialize
    local idx_s::Array{Int64}     # initialize



    # Compute the total s (current position along track)
    s_total = s % s_target

    # Compute the index
    DistS = ( s_total - oldS ).^2

    idx_s = findmin(DistS,1)[2]              # contains both indices for the closest distances for both oldS !!

    off = 0
    idx_s = idx_s + off
    
    
    selectedStates.selStates[i = 1:Np,j=1:4] = oldTraj.oldTraj[j=idx_s[1]:1:idx_s[1]+Np-1,i=1:4,selected_laps[1]]
    selectedStates.selStates[i=Np+1:2*Np,i=1:4] = oldTraj.oldTraj[i=idx_s[2]-N_points:1:idx_s[2]+Np-N_points-1,i=1:4,selected_laps[2]]  
    

    selectedStates.statesCost[i=1:Np] = oldTraj.cost2Target[j=idx_s[1]:idx_s[1]+Np-1] 
    selectedStates.statesCost[i=Np+1:2*Np] = oldTraj.cost2Target[j=idx_s[2]-N_points:idx_s[2]-N_points+Np-1]

    println("current s = ",s)
    #println("selected states = ",selectedStates.selStates[:,6])
    #println("states' cost = ",selectedStates.statesCost)
    #println("check traj = ", oldTraj.oldTraj[idx_s[1]:idx_s[1]+12,:,selected_laps[1]])


    nothing
end
