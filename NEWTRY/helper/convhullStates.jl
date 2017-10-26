# This function selects, from the previous laps, the states that will be used to build the convex hull in the LMPC formulation. 
# Moreover, it returns their costs, as evaluated in the previous laps. 

# z[1] = s 
# z[2] = ey
# z[3] = epsi
# z[4] = v


function convhullStates(oldTraj::classes.OldTrajectory, posInfo::classes.PosInfo, mpcParams::classes.MpcParams,lapStatus::classes.LapStatus,selectedStates::classes.SelectedStates, obs::Array{Float64},modelParams::classes.ModelParams,obstacle::classes.Obstacle,simVariables::classes.SimulationVariables)

    # Read Inputs
    s               = posInfo.s          # current s 
    s_target        = posInfo.s_target   # s value of the target

    # Parameters
    Np              = selectedStates.Np  # number of states to choose per each last lap
    Nl              = selectedStates.Nl  # number of previous laps to take in consideration for the convex hull
    N               = mpcParams.N        # prediction horizon    
    dt              = modelParams.dt     # time step
    n_pf            = simVariables.n_pf  # number of path following laps 
    r_s             = obstacle.r_s       # radius on the s coordinate of the ellipse describing the obstacle
    r_ey            = obstacle.r_ey      # radius on the ey coordinate of the ellipse describing the obstacle
   
    
 
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

    # Propagate the obstacle for the prediction horizon

    obs_prop_s  = obs[1,1,:] + dt*N*obs[1,3,:]
    obs_prop_ey = obs[1,2,:]
   

    
    
    for j = 0:(Nl-1)

        selectedStates.selStates[i=(j*Np)+1:(j+1)*Np,i=1:4] = oldTraj.oldTraj[i=idx_s[j+1]-(j*N_points):idx_s[j+1]+Np-(j*N_points)-1,i=1:4,selected_laps[j+1]]  # select the states from lap j...
        
        selectedStates.statesCost[i=(j*Np)+1:(j+1)*Np] = oldTraj.cost2target[i=idx_s[j+1]-(j*N_points):idx_s[j+1]-(j*N_points)+Np-1,selected_laps[j+1]]  # and their cost

        if obstacle.lap_active == true   # if the obstacles are on the track, check if any of the selected states interferes with the propagated obstacle

            for n=1:obstacle.n_obs
                ellipse_check = (((selectedStates.selStates[i=(j*Np)+1:(j+1)*Np,1]-obs_prop_s[n])/r_s)^2) + (((selectedStates.selStates[i=(j*Np)+1:(j+1)*Np,2]-obs_prop_ey[n])/r_ey)^2)
                
                if any(x->x<=1, ellipse_check) == true  # if any of the selected states is in the ellipse

                    index = find(ellipse_check.<=1)     # find all the states in the ellipse 

                    mpcParams.Q_obs[i=(j*Np)+(index[1]-obstacle.inv_step)+1:(j+1)*Np] =  10   # and set the values of the weight to 10, so that they are excluded from optimization
                end
            end
        end     
        
    end



    nothing
end
