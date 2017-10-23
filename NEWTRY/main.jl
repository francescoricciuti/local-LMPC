# This is the main file, which should be launched to execute the program. To launch this file, navigate to your local directory in 
# which this file is saved, execute Julia, include this file ( typing 'include("main.jl")'), and run it

using JuMP
using Ipopt
using JLD

#### Include Helpers
include("helper/classes.jl")
include("helper/initializeParameters.jl")
include("helper/initializeModels.jl")
include("helper/prova.jl")
include("helper/solveMPCmodels.jl")
include("helper/simModels.jl")
include("helper/saveOldTraj.jl")
include("helper/convhullStates.jl")
include("helper/obstacle.jl")




#### Load Maps
println("loadMap.......")
include("helper/loadTestMap.jl")
x_track, y_track = loadTestMap() #load a track to test controller
println("loaded")

#### Define Types (from helper/classes.jl)
println("define types ........")

lapStatus         = classes.LapStatus()
oldTraj           = classes.OldTrajectory()
selectedStates    = classes.SelectedStates()
mpcParams         = classes.MpcParams()
posInfo           = classes.PosInfo()
mpcSol            = classes.MpcSol()
trackCoeff        = classes.TrackCoeff()
modelParams       = classes.ModelParams()
simVariables      = classes.SimulationVariables()
obstacle          = classes.Obstacle()


#### Initialize states 

z_Init    = zeros(6)
z_Init[1] = 0                  #x
z_Init[2] = 0                  #y
z_Init[3] = 0.4*cos(z_Init[5]) #v_x
z_Init[4] = 0.4*sin(z_Init[5]) #v_y
z_Init[5] = 0                  #psi
z_Init[6] = 0                  #psi_dot

#### Initialize parameters

InitializeParameters(mpcParams,trackCoeff,modelParams,oldTraj,mpcSol,lapStatus,simVariables,selectedStates,obstacle)

#### Calculate position of s_target

@show posInfo.s_target = (size(x_track)[2]-1)*trackCoeff.ds

#### Define simulation parameters
dt         = modelParams.dt                   # time step
buffersize = simVariables.buffersize          # buffersize ( to have consistent dimensions in all the arrays of the code)
n_laps     = simVariables.n_laps              # number of laps we want to simulate
t          = collect(0:dt:(buffersize-1)*dt)  # time array
postbuff   = simVariables.postbuff            # number of postbuff iteration to save after end of lap
s_target   = posInfo.s_target                 # position of the target
dynModel   = simVariables.dynModel            # boolean variable telling the simulator which model to use
Nl         = selectedStates.Nl                # number of previous laps to consider in the convex hull
n_obs      = obstacle.n_obs                   # number of obstacles in the track
obs_detect = obstacle.obs_detect              # maximum distance at which we can detect obstacles

#### Initialize Models

println("Initialize Models........")
mdl_Path = initPathFollowingModel(mpcParams,modelParams,trackCoeff)
mdl_LMPC = initLearningModel(mpcParams,modelParams,trackCoeff,selectedStates)
mdl_obs_LMPC=initObsModel(mpcParams,modelParams,trackCoeff,selectedStates,obstacle)



#### Create log variables (needed for later plotting)

selStates_log  = zeros(Nl*selectedStates.Np,4,length(t),n_laps)   #array to log the selected states in every iteration of every lap
statesCost_log = zeros(Nl*selectedStates.Np,length(t),n_laps)     #array to log the selected states' costs in every iteration of every lap
z_pred_log     = zeros(mpcParams.N+1,4,length(t),n_laps)          #array to log the predicted states for all time steps for later plotting
u_pred_log     = zeros(mpcParams.N,2,length(t),n_laps)            #array to log the predicted inputs for all time steps for later plotting
cost_log       = zeros(7,length(t),n_laps)                        #logs the MPC costs with values for all expresions in the cost function (terminal cost, lane cost, input cost etc.)
alpha_log      = zeros(Nl*selectedStates.Np,length(t),n_laps)     #logs the coefficients for the convex hull
obs_log        = zeros(length(t),3,n_obs,n_laps)                  #logs the info about the obstacle 

#### Initialize variables needed for the main loops 

j              = 1                                               # set j to one as we start at the first lap 
z_final_x      = zeros(1,4)::Array{Float64,2}                    # initialize final states...
u_final        = zeros(1,2)::Array{Float64,2}                    # ...and final control inputs. At the end of every lap are used to initialize the car in the new lap
obs_final      = zeros(1,3,n_obs)::Array{Float64,3}              # save last obstacles states of one laps to initialize the successive one
tt             = zeros(length(t),1)::Array{Float64,2}            # logs the time used by the solveLearningMpcProblem
tt_it          = zeros(length(t),1)::Array{Float64,2}            # logs the total compuational time for one iteration of the inner loop
zCurr_s        = zeros(length(t),4)::Array{Float64,2}            # (s, ey, epsi, v) every time an iteration ends, save the current states here so that, at the end of a given lap, we can save them in OldTrajectory
zCurr_x        = zeros(length(t),6)::Array{Float64,2}            # (x, y, v_x,v_y, psi, psi_dot) same as above, in x-y frame
uCurr          = zeros(length(t),2)::Array{Float64,2}            # every time an iteration ends, save here the applied input so that, at the end of a given lap, we can save them
curvature_curr = zeros(length(t))::Array{Float64,1}              # curvature as calculated in every iteration of the inner loop
obs_curr       = zeros(length(t),3,n_obs)::Array{Float64,3}      # info about the obstacle in the current lap

#### Set the initial conditions as zCurr_x to initialize first iteration of first lap

zCurr_x[1,1]   = z_Init[1]  
zCurr_x[1,2]   = z_Init[2] 
zCurr_x[1,3]   = z_Init[3]
zCurr_x[1,4]   = z_Init[4]
zCurr_x[1,5]   = z_Init[5]
zCurr_x[1,6]   = z_Init[6]

#### Set initial conditions on the obstacles

obs_curr[1,1,:] = obstacle.s_obs_init
obs_curr[1,2,:] = obstacle.ey_obs_init
obs_curr[1,3,:] = obstacle.v_obs_init


#### Start the main loop 

for j=1:n_laps                 # main loop over all laps

    

    if j == obstacle.lap_active            # if its time to put the obstacles in the track
        obstacle.obstacle_active = true    # tell the system to put the obstacles on the track
    end

    if j > obstacle.lap_active             # initialize current obstacle states with final states from the previous lap
        obs_curr[1,:,:] = obs_final

    end
    

    no_solution_found = 0      #counts number of unsuccesful attempts. If too many unsuccesful attempts, the lap will be terminated

    lapStatus.currentLap = j   # keep track of the current lap updating the value in lapStatus

    # if we are not in the first lap, initialize states and inputs with the last values obtained from the previous lap
    if j>1    
        zCurr_x = zeros(length(t),6)
        zCurr_x[1,:] = z_final_x
        uCurr[1,:] = u_final
    end


#### Start the inner loop, i.e. iterations in each lap
    
    finished            = false     # it is set to true if finish line is reached and round is finished
    i                   = 1         # index for iterations in one round
    lapStatus.currentIt = 1         # since we are starting a new lap, current iteration is 1

   while i<=length(t)-1 && !finished    # as long as we have not reached the maximal iteration time for one round or ended the round
   #for indice=1:10
        tic()  # tic for iteration time calculation (tt_it)
        mpcParams.Q_obs = ones(Nl*selectedStates.Np)

        if j > 1
            if i == (postbuff+2)
                oldTraj.oldTraj[oldTraj.costLap[j-1]+1:oldTraj.costLap[j-1]+postbuff+1,2:4,j-1] = zCurr_s[1:postbuff+1,2:4]
                oldTraj.oldTraj[oldTraj.costLap[j-1]+1:oldTraj.costLap[j-1]+postbuff+1,1,j-1] = zCurr_s[1:postbuff+1,1] + s_target
            end
        end



        # transforms states form x-y coordinates to s-ey coordinates and approximates a curvature around the current postion, returning corresponding polynomial coefficients
        zCurr_s[i,:], trackCoeff.coeffCurvature = trackFrameConversion(zCurr_x[i,:],x_track,y_track,trackCoeff, i,oldTraj,j,true)

        if i == 1 && zCurr_s[1,1] > 2 # if we are in the first iteration and, after the conversion from x-y to s-ey, it results that the s coordinate is greater than 2 
                                      #(meaning, for example, that we haven't crossed the start line yet, thus having an s coordinate big), we force the s coordinate to
                                      #zero (thus forcing it to be on the start line)
            zCurr_s[1,1]= 0
        end

        posInfo.s         = zCurr_s[i,1]    # save the s coordinate of the current position. 
        curvature_curr[i] = trackCoeff.coeffCurvature[1]*posInfo.s^3+trackCoeff.coeffCurvature[2]*posInfo.s^2+trackCoeff.coeffCurvature[3]*posInfo.s+trackCoeff.coeffCurvature[4] # calculate the curvature just for debugging
        
        #### Check if we have reached the finish line
        
        if zCurr_s[i,1] >= posInfo.s_target #if the car has crossed the finish line
            println("Reached finish line at step $i")
            finished = true                 # tell the system that we finished current lap. "finished=true" stops the while loop
            
            break

        elseif i >1 && zCurr_s[i-1,1] >= posInfo.s_target-1 && zCurr_s[i,1]<1 # on a closed track the finish line lies close to the starting point. It can happen that the discrete calculated s postion is bigger than the end of the track and already in the new round. 
                                                                              # To still detect the finish line this statement is added
            println("alternative detection of finish line")
            println("Reached finish line at step $i")
            finished = true                 # tell the system that we finished current lap. "finished=true" stops the while loop
            
            break
        end   

        #### If the obstacles are on the track, find the nearest one

        if obstacle.obstacle_active == true

            # if posInfo.s_target-posInfo.s < obs_detect  # meaning that I could possibly detect obstacles after the finish line

            #     index1=find(obs_curr[i,1,:].< obs_detect+posInfo.s-posInfo.s_target)  # look for obstacles that could cause problems

            dist,index=findmin(obs_curr[i,1,:]-posInfo.s)     # find the minimum distance and the index of the nearest obstacle_active
            println("dist= ",dist)
            obs_near = obs_curr[i,:,index]
        end


        #### Solve the MPC problem


        n_pf = simVariables.n_pf  # number of path following laps

        tic() # tic for learning MPC time calculations (tt)
        if j <= n_pf       # if we have not completed the path following laps yet, just do path following

            solvePF_MPC(mdl_Path,mpcSol,mpcParams,trackCoeff,modelParams,zCurr_s[i,:]',uCurr[i,:]')

###########################################################################################################################################
        elseif j > n_pf    # if we have already completed all the path following laps, compute the states needed for the convex hull and solve the LMPC

            convhullStates(oldTraj, posInfo, mpcParams,lapStatus, selectedStates, obs_curr[i,:,:],modelParams,obstacle,simVariables)

            if obstacle.obstacle_active == false 
                solveLearning_MPC(mdl_LMPC,mpcSol,mpcParams,trackCoeff,modelParams,zCurr_s[i,:]',uCurr[i,:]',selectedStates)
            elseif obstacle.obstacle_active == true
                if dist > obs_detect || zCurr_s[i,1]>=obs_near[1,1]+obstacle.r_s
                    solveLearning_MPC(mdl_LMPC,mpcSol,mpcParams,trackCoeff,modelParams,zCurr_s[i,:]',uCurr[i,:]',selectedStates)
                    println("FLAG LMPC")
                    # println("current s state at iteration $i= ",zCurr_s[i,:])
                    # println("current xy state at it $i= ",zCurr_x[i,:])
                    # println("applied control at it $i= ",[mpcSol.a_x mpcSol.d_f])
                    # println("predicted trajectory at it $i= ",mpcSol.z)
                elseif dist <= obs_detect && zCurr_s[i,1]<obs_near[1,1]+obstacle.r_s
                    solveObs_LMPC(mdl_obs_LMPC,mpcSol,mpcParams,trackCoeff,modelParams,zCurr_s[i,:]',uCurr[i,:]',selectedStates,obs_near,obstacle)
                    println("FLAG OBS")
                   
                end
            end

            alpha_log[:,i,j] = mpcSol.alpha # save the coefficients for convex hull computed in iteration i of lap j
        
###########################################################################################################################################à
        end
        tt[i]= toq() # toq for learning MPC time calculation 
 
        uCurr[i,:]      = [mpcSol.a_x mpcSol.d_f] # set as current input the first of the inputs computed bt the MPC and apply it to the simulator ("simModel")

        if dynModel == true
            zCurr_x[i+1,:]  = simModel_exact_dyn_x(zCurr_x[i,:],uCurr[i,:],modelParams)  # simulates the dynamical model with the inputs generated by the controller
        else 
            zCurr_x[i+1,:]  = simModel_kin_x(zCurr_x[i,:],uCurr[i,:],modelParams) # simulates the kinematic model with the inputs generated by the controller
        end

        if obstacle.obstacle_active == true
            obs_curr[i+1,:,:] = obstaclePosition(obs_curr[i,:,:],modelParams,obstacle)
        end

        #### Once the MPC is solved, save the solution in the log variables for later plotting
       
        cost_log[:,i,j]               = mpcSol.cost  # save the optimal cost resulting from iteration i of lap j
        z_pred_log[:,:,i,j]           = mpcSol.z     # save the states computed in iteration i of lap j
        u_pred_log[:,:,i,j]           = mpcSol.u     # save the control actions computed in iteration i of lap j
        selStates_log[:,:,i,j]        = selectedStates.selStates 
        statesCost_log[:,i,j]         = selectedStates.statesCost
        obs_log[i,:,:,j]              = obs_curr[i,:,:]
        
        tt_it[i]= toq() # toq for iteration time calculation 

        #### Print on the screen the solution of the current optimization iteration
        println(" Lap $j , Solving step $i of $(length(t)), s = $(zCurr_s[i,1]) - Status: $(mpcSol.solverStatus)")

        i = i + 1                 # count up i for next iteration in this round
        lapStatus.currentIt = i   # update the value of the current iteration

        println("Current it: ",i) # print on the screen the current iteration
    end

    #### When here, either the finish line has been crossed, or we have reached the maximum number of iterations allowed (i.e. buffersize)

    # Comunicate if we exceeded the maximum number of allowed iterations
    if i >= length(t)
        println("used whole length t. Either adapt array size to track or problems in solving process")
    end

    z_final_x = zCurr_x[i,:]  # this is the state in which we arrived by applying the last control input computed by the MPC
    u_final   = uCurr[i-1,:]  # this is the last control input computed by the MPC. Note that we have to index with i-1 because 
                              # in the meanwhile i has been increased by 1. In case this passage is not clear, just go where simModel_exact_dyn_x
                              # was solved and everything will be clear.
    obs_final = obs_curr[i,:,:]
    println("=================\nFinished Solving. Avg. time = $(mean(tt[1:i])) s")
    println("Finished Lap Nr. $j")

    i = i-1                     # we check if we have reached the finish line at the start of the inner loop, but we increased the iteration at the end of
                                # the previous inner loop. This means that when we realize that we have reached the finish line and we break the while loop
                                # we already increased the iteration number of one, thus saving in currentIt an iteration that we have never done.
                                # For this reason we have to decrease of one the iteration number...
    lapStatus.currentIt = i     # ... and also the currentIt value. This way currentIt will be the number of iterations we have done to reach the finish line

    #### Save data of the current lap in OldTrajectory

    saveOldTraj(oldTraj,cost_log,zCurr_s, zCurr_x,uCurr,lapStatus,simVariables,z_pred_log,u_pred_log,modelParams,curvature_curr,alpha_log,posInfo)

end

#set filename to save data
#### numbering to generate results to keep 
filename = string("data/"string(Dates.today()),"-",Dates.format(now(), "HH-MM"),"-Data.jld")
if isfile(filename)
    filename = string("data/"string(Dates.today()),"-",Dates.format(now(), "HH-MM"),"-Data-2.jld")
    warn("File already exists. Added extension \"-2\" ")
end
println("Save data to $filename .......")


#safe data to file
jldopen(filename, "w") do file
    #addrequire(file, classes) #ensures that custom data types are working when loaded

    JLD.write(file, "x_track", x_track)
    JLD.write(file, "y_track", y_track)
    JLD.write(file, "trackCoeff", trackCoeff)
    JLD.write(file, "modelParams", modelParams)
    JLD.write(file, "mpcParams", mpcParams)
    JLD.write(file, "buffersize", buffersize)
    JLD.write(file, "oldTraj", oldTraj)
    JLD.write(file, "selectedStates", selStates_log)
    JLD.write(file, "obstacleData", obstacle)
    JLD.write(file, "obs", obs_log)
    

end
println("finished")

