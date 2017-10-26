# This file contains all the functions related to the obstacles

# this function computes the position of the obstacle
function obstaclePosition(obs_now::Array{Float64},modelParams::classes.ModelParams,obstacle::classes.Obstacle)

    n_obs = obstacle.n_obs
    obs_next=zeros(1,3,n_obs)::Array{Float64}    # initialize the array that will contain the next states of the obstacle
    
    dt = modelParams.dt

    obs_next[:,1,:] = obs_now[:,1,:] + obs_now[:,3,:]*dt # s
    obs_next[:,2,:] = obs_now[:,2,:]                     # ey
    obs_next[:,3,:] = obs_now[:,3,:]                     # v

    return obs_next

end  

