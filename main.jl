using JuMP
using Ipopt
using JLD

#### Include Helpers
include("helper/classes.jl")
include("helper/functions.jl")
include("helper/initializeModel.jl")
include("helper/convhullStates.jl")
include("helper/solveMpcProblem.jl")
include("helper/simModel.jl")
include("helper/localizeVehicleCurvAbs.jl")
include("helper/computeObstaclePos.jl")
include("helper/adjustSafeSet.jl")

#### Load Maps 

#just loads one specified track with distances betwwen points =1 m and returnx x and y coordinatates
println("loadMap.......")
include("helper/loadTestMap.jl")
x_track, y_track = loadTestMap() #load a track to test controller
println("loaded")

#### Define the Types

println("define types ........")
oldTraj                     = classes.OldTrajectory()
posInfo                     = classes.PosInfo()
lapStatus                   = classes.LapStatus(1,1)
mpcSol                      = classes.MpcSol()
trackCoeff                  = classes.TrackCoeff()      # info about track (at current position, approximated)
modelParams                 = classes.ModelParams()
mpcParams                   = classes.MpcParams()
selectedStates              = classes.SelectedStates()
mpcCoeff                    = classes.MpcCoeff()

buffersize                  = 1101 # buffersize has to be big enough to fit all iterations for the path following trajectory
n_rounds                    = 3   # how many rounds we want to simulate


#### Initialize States

z_Init    = zeros(6)
z_Init[1] = 0                  #x
z_Init[2] = 0                  #y
z_Init[3] = 0.4*cos(z_Init[5]) #v_x
z_Init[4] = 0.4*sin(z_Init[5]) #v_y
z_Init[5] = 0                  #psi
z_Init[6] = 0.0                #psi_dot


#### Initialize Parameters

InitializeParameters(mpcParams,trackCoeff,modelParams,oldTraj,mpcCoeff,mpcSol,lapStatus,buffersize,selectedStates) #initillize value for the classes specified above

#### Calculate position of s_target

@show posInfo.s_target = (size(x_track)[2]-1)*trackCoeff.ds

#### Initialize Model

println("Initialize Model........")
#the controller relies on two independent models for the pathfollowing rounds and all following learning mpc rounds
mdl_Path = initPathFollowingModel(mpcParams,modelParams,trackCoeff,mpcCoeff,oldTraj.n_oldTraj)
mdl_LMPC = initLearningModel(mpcParams,modelParams,trackCoeff,mpcCoeff,oldTraj.n_oldTraj,selectedStates)

println("Initial solve........")
solve(mdl_LMPC.mdl)# intial solve is necessary for LMPC model to prevent an invalid number error. seems to be caused by problems with the coefficients of the terminal cost and terminal set in the controller.
println("Initial solve done!")
println("*******************************************************")

# Simulation parameters
dt                          = modelParams.dt
t                           = collect(0:dt:(buffersize-1)*dt)

Pcurvature = zeros(length(t),2) # just for debugging process. Collects the approximated curvature at all positions from the localizing function. 


trackCoeff.coeffCurvature   = [0.0;0.0;0.0;0.0;0.0]  # polynomial coefficients for curvature approximation (zeros for straight line)
trackCoeff.nPolyCurvature = 4                        # has to be 4, cannot be changed freely. At the moment orders are still hardcoded in some parts of localizeVehicleCurvAbs
trackCoeff.nPolyXY = 6                               # has to be 6, cannot be changed freely. At the moment orders are still hardcoded in some parts of localizeVehicleCurvAbs

# create log variables for later plotting

selStates_log  = zeros(2*selectedStates.Np,4,length(t),n_rounds)  #array to log the selected states in every iteration of every round
statesCost_log = zeros(2*selectedStates.Np,length(t),n_rounds)    #array to log the selected states' costs in every iteration of every round
z_pred_log     = zeros(mpcParams.N+1,4,length(t),n_rounds)        #array to log the predicted states for all time steps for later plotting
u_pred_log     = zeros(mpcParams.N,2,length(t),n_rounds)          #array to log the predicted inputs for all time steps for later plotting
#lambda_log     = zeros(oldTraj.n_oldTraj,length(t),n_rounds)      #array to log the predicted lambdas
cost           = zeros(5,length(t),n_rounds)                      #logs the MPC costs with values for all expresions in the cost function (terminal cost, lane cost, input cost etc.)


j = 1                                     # set j to one as we start at the first round 
z_final_x = zeros(1,4)::Array{Float64,2}  # initialize final states...
u_final = zeros(1,2)                      # ...and final control inputs. At the end of every round are used to initialize the car in the new round

#### Start the main loop 

for j=1:n_rounds     #loop over all rounds
    
    no_solution_found = 0 #counts number of unsuccesful attempts. If too many unsuccesful attempts, the round will be terminated

    lapStatus.currentLap = j
    tt          = zeros(length(t),1)     #logs the time used by the solveLearningMpcProblem
    #tt1         = zeros(length(t),1)     #logs the time used by safe set operations ( addTrajectory, deleteTrajectory, interpolateTrajectory)
    tt_total    = zeros(length(t),1)     #logs the total compuational time for one iteration of the inner loop
    zCurr_s     = zeros(length(t),4)          # s, ey, epsi, v
    zCurr_x     = zeros(length(t),6)          # x, y, v_x,v_y, psi, psi_dot
    uCurr       = zeros(length(t),2)
    curvature_curr = zeros(length(t))

    #setup point for vehicle on track in first round. gets overwritten in other rounds

    zCurr_x[1,1] = z_Init[1] 
    zCurr_x[1,2] = z_Init[2] 
    zCurr_x[1,3] = z_Init[3]
    zCurr_x[1,4] = z_Init[4]
    zCurr_x[1,5] = z_Init[5]
    zCurr_x[1,6] = z_Init[6]

    # if we are not in the first round, initialize states and imputs with the values obtained from the previous optimization
    if j>1       
        zCurr_x[1,:] = z_final_x
        uCurr[1,:] = u_final
    end

#### Start the inner loop, i.e. iterations in each round
    
    finished = false     # it is set to true if finish line is reached and round is finished
    i        = 1         # index for iterations in one round

    while i<=length(t)-1 && !finished    # as long as we have not reached the maximal iteration time for one round or ended the round

        tic()  # tic for overall time calculation (tt_total)

        # transforms states form x-y coordinates to s-ey coordinates and approximates a curvature around the current postion, returning corresponding polynomial coefficients
        zCurr_s[i,:], trackCoeff.coeffCurvature = localizeVehicleCurvAbs(zCurr_x[i,:],x_track,y_track,trackCoeff, i, mpcParams.N, modelParams.dt, Pcurvature)

        if i == 1 && zCurr_s[1,1] > 2 # if the xy coordinates from the loaded states are closest to a point at the end of the track not at the beginning of the track in the first step the s coordinate is reset to zero
            # warn("closest point was before finish line, forced s =0")
            # println("x:$(zCurr_x[i,1]), y:$(zCurr_x[i,2])")
            zCurr_s[1,1]= 0
        end

        posInfo.s   = zCurr_s[i,1]
        curvature_curr[i] = trackCoeff.coeffCurvature[1]*posInfo.s^4+trackCoeff.coeffCurvature[2]*posInfo.s^3+trackCoeff.coeffCurvature[3]*posInfo.s^2+trackCoeff.coeffCurvature[4]*posInfo.s +trackCoeff.coeffCurvature[5] # calculate the curvature just for debugging
        
        if zCurr_s[i,1] >= posInfo.s_target #if the car has crossed the finish line
            println("Reached finish line at step $i")
            finished = true
            #we count up here as the first round ends just as we would do if the loop get terminiated because i is >= length(t). we count it down later on again to get right index
            i = i + 1
            lapStatus.currentIt = i
            break
        elseif i >1 && zCurr_s[i-1,1] >= posInfo.s_target-1 && zCurr_s[i,1]<1 # on a closed track the finish line lies close to the starting point. it can be that the discrete calculated s postion is bigger than the end of the track and already in the new round. to still detect the finish line we added this statement
            println("alternative detection of finish line")
            println("Reached finish line at step $i")
            finished = true
            #we count up here as the first round ends just as we would do if the loop get terminiated because i is >= length(t). we count it down later on again to get right index
            i = i + 1
            lapStatus.currentIt = i
            break
        end         

        if j > 3

            convhullStates(oldTraj, mpcCoeff, posInfo, mpcParams,lapStatus, selectedStates)

        end


        tic() # tic for learning MPC time calculations (tt)

        if j < 4 # use path following if we are in the first 3 rounds
            solvePathFollowMpc!(mdl_Path,mpcSol,mpcCoeff,mpcParams,trackCoeff,lapStatus,posInfo,modelParams,zCurr_s[i,:],[mpcSol.a_x;mpcSol.d_f],i)
        elseif j >= 4 
            solstat = solveLearningMpcProblem!(mdl_LMPC,mpcSol,mpcCoeff,mpcParams,trackCoeff,lapStatus,posInfo,modelParams,zCurr_s[i,:],[mpcSol.a_x;mpcSol.d_f],i,selectedStates)
            if solstat == false
                no_solution_found += 1
            else
                no_solution_found = 0
            end
            if no_solution_found >= 3
                warn("Over 3 unsuccesful iterations. Abort solving!")
                break
            end    
            
        end
        tt[i]= toq()

        uCurr[i,:]  = [mpcSol.a_x mpcSol.d_f] 
        zCurr_x[i+1,:]  = simModel_exact_dyn_x(zCurr_x[i,:],uCurr[i,:],modelParams.dt,modelParams)  # simulates the model with the inputs generated by the controller

        # save the solution of the controller for future plotting
       
        cost[:,i,j]               = mpcSol.cost
        #lambda_log[:,i,j]         = mpcSol.lambda
        z_pred_log[:,:,i,j]       = mpcSol.z
        u_pred_log[:,:,i,j]       = mpcSol.u
        #ssInfOn_log[:,i,j]       = mpcSol.ssInfOn
        selStates_log[:,:,i,j]    = selectedStates.selStates 
        statesCost_log[:,i,j]     = selectedStates.statesCost

        tt_total[i]= toq()


         if string(mpcSol.solverStatus) != "Optimal" #print every time the solver fails to find an optimal solution
            println(" Time: $(tt[i]) s, Solving step $i of $(length(t)), s = $(zCurr_s[i,1]) - Status: $(mpcSol.solverStatus)")
        end
        

        i = i + 1#count up i for next iteration in this round
        lapStatus.currentIt = i

        println("Current it: ",i)
    end
    #at this point either the finish line is crossed or the array to safe data is full. 
    if i >= length(t)
        println("used whole length t. either adapt array size to track or problems in solving process")
    end

    i = i-1 
    lapStatus.currentIt -= 1 # has finished lap already so we need to count both coutners one down as we counted them up after we have crossed the finish line
    z_final_x = zCurr_x[i,:]
    u_final = uCurr[i,:]
    println("=================\nFinished Solving. Avg. time = $(mean(tt[1:i])) s")
    println("Finished Lap Nr. $j")

    # Save states in oldTraj:
    # --------------------------------
    tic()

    #saves ego vehicle data into oldTraj for the safe set and for future plotting
    saveOldTraj(oldTraj,zCurr_s, zCurr_x,uCurr,lapStatus,buffersize,modelParams.dt, cost[:,:,j],z_pred_log[:,:,:,j],u_pred_log[:,:,:,j], mpcSol, curvature_curr)
    tt3= toq()

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
    #JLD.write(file, "obstacle", obstacle)
    JLD.write(file, "modelParams", modelParams)
    JLD.write(file, "mpcParams", mpcParams)
    JLD.write(file, "buffersize", buffersize)
    JLD.write(file, "oldTraj", oldTraj)
    JLD.write(file, "mpcCoeff",mpcCoeff)

end
println("finished")
