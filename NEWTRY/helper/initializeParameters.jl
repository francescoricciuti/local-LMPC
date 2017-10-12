# This file initializes all the parameters 

function InitializeParameters(mpcParams::classes.MpcParams,trackCoeff::classes.TrackCoeff,modelParams::classes.ModelParams,
                                oldTraj::classes.OldTrajectory,mpcSol::classes.MpcSol,lapStatus::classes.LapStatus,simVariables::classes.SimulationVariables,selectedStates::classes.SelectedStates)
    
    simVariables.buffersize     = 2000      # used to initialize the dimensions of the variables in which we will save the data of the Simulations 
    buffersize                  = simVariables.buffersize
    simVariables.n_laps         = 10       # number of laps we want to simulate 
    simVariables.n_pf           = 3        # number of path following laps (must be at least 2)
    simVariables.postbuff       = 30       # number of postbuffer iteration to save
    dynModel                    = false    # boolean variable to tell the simulator which model to use (dynModel=True-->it'll use dynamic model, dynModel=False-->it'll use kinematic model)

    mpcParams.N                 = 17                        #lenght of prediction horizon
    mpcParams.vPathFollowing    = 0.6                       # reference velocity for the path following stage
    mpcParams.QderivZ           = 0.1*[0,0.0,0.1,0.1]       # weights for the states in the derivative cost
    mpcParams.QderivU           = 0.1*[1.0,10]              # weights for the control inputs in the derivative cost
    mpcParams.R                 = 0.0*[1.0,1.0]             # weights on the control inputs
    mpcParams.Q                 = [0.0,50.0,0.1,10.0]       # weights on the states for path following
    mpcParams.Q_cost            = 0.7                       # weight on the cost to get from a given point to the targe
    mpcParams.Q_lane            = 0.9                       # weight on the soft constraint on the Q_lane
    mpcParams.Q_alpha           = 1.0                       # weight on the soft constraint for convex hull

    trackCoeff.nPolyCurvature   = 4                       # 4th order polynomial for curvature approximation
    trackCoeff.nPolyXY          = 6 
    trackCoeff.coeffCurvature   = [0.0;0.0;0.0;0.0;0.0]         # polynomial coefficients for curvature approximation (zeros for straight line)
    trackCoeff.width            = 2                   # width of the track (0.6m)
    trackCoeff.ds               = 1//10#4//100#1//10 # is defined as a rational number so we can use it to calculate indices in matrix. with float becomes error

    modelParams.l_A             = 0.125
    modelParams.l_B             = 0.125 
    modelParams.dt              = 0.1
    modelParams.u_lb            = ones(mpcParams.N,1) * [-0.6  -pi/4]  #-0.6 for braking    # lower bounds on control actions
    modelParams.u_ub            = ones(mpcParams.N,1) * [ 1.2   pi/4]       #1.2           # upper bounds on control actions
    modelParams.z_lb            = ones(mpcParams.N+1,1) * [-Inf -Inf -Inf  0]           # lower bounds on states
    modelParams.z_ub            = ones(mpcParams.N+1,1) * [ Inf  Inf  Inf  1.5]            # upper bounds on states
    modelParams.v_max           = 2.0                                                      # maxium velocity for trajectory tracking
    modelParams.max_alpha       = 10                                                       # maximum slip angle 
    modelParams.mass            = 1.98 # kg
    modelParams.mu              = 0.85
    modelParams.g               = 9.81 # m/s^2
    modelParams.I_z             = 0.03 # kg * m^2
    modelParams.B               = 6.0  
    modelParams.C               = 1.6  

    selectedStates.Np           = 20                            # Number of points to take from each previous trajectory to build the convex hull
    selectedStates.selStates    = zeros(2*selectedStates.Np,4)  
    selectedStates.statesCost   = zeros(2*selectedStates.Np)

    oldTraj.n_oldTraj           = simVariables.n_laps                                                     # number of old Trajectories for safe set
    oldTraj.oldTraj             = zeros(buffersize,4,oldTraj.n_oldTraj)                  # old trajectories in s-ey frame
    oldTraj.oldTrajXY           = zeros(buffersize,6,oldTraj.n_oldTraj)                  # old trajectories in x-y frame
    oldTraj.oldNIter            = zeros(oldTraj.n_oldTraj)                               # how many iterations to complete that particular iteration
    oldTraj.oldInput            = zeros(buffersize,2,oldTraj.n_oldTraj)                  # old Inputs 
    oldTraj.costs               = zeros(5,buffersize,oldTraj.n_oldTraj)                  # optimal values of the cost function elements for each optimization iteration
    oldTraj.z_pred_sol          = zeros(mpcParams.N+1,4,buffersize,oldTraj.n_oldTraj)    # predicted states (in s-ey frame) for each iteration of past rounds
    oldTraj.u_pred_sol          = zeros(mpcParams.N,2,buffersize,oldTraj.n_oldTraj)      # predicted input for each iteration of past rounds
    oldTraj.cost2target         = zeros(buffersize,oldTraj.n_oldTraj)                    # number of iterations needed to arrive at the target
    oldTraj.curvature           = zeros(buffersize,oldTraj.n_oldTraj)                    # all the curvatures calculated in each iteration of each lap
    oldTraj.oldAlpha            = zeros(2*selectedStates.Np,buffersize,oldTraj.n_oldTraj)# all the alphas from all iterations of all LMPC laps
    oldTraj.costLap             = zeros(oldTraj.n_oldTraj)


    mpcSol.u                    = zeros(mpcParams.N,2)              # array containing all the control inputs computed by the MPC at a given iteration
    mpcSol.z                    = zeros(mpcParams.N+1,4)            # array containing all the states computed by the MPC at a given iteration
    mpcSol.cost                 = zeros(5)                          # optimal costs as computed by the MPC at a given iteration
    mpcSol.alpha                = zeros(buffersize,2*selectedStates.Np,oldTraj.n_oldTraj)  # coefficients of the convex hull
    
    lapStatus.currentLap        = 1         # initialize lap number
    lapStatus.currentIt         = 0         # current iteration in lap 

    

    
end