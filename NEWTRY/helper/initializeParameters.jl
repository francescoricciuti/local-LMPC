# This file initializes all the parameters 

function InitializeParameters(mpcParams::classes.MpcParams,trackCoeff::classes.TrackCoeff,modelParams::classes.ModelParams,
                                oldTraj::classes.OldTrajectory,mpcSol::classes.MpcSol,lapStatus::classes.LapStatus,simVariables::classes.SimulationVariables,selectedStates::classes.SelectedStates,obstacle::classes.Obstacle)
    
    simVariables.buffersize     = 2000      # used to initialize the dimensions of the variables in which we will save the data of the Simulations 
    buffersize                  = simVariables.buffersize
    simVariables.n_laps         = 9       # number of laps we want to simulate 
    simVariables.n_pf           = 5        # number of path following laps (must be at least Nl)
    simVariables.postbuff       = 40       # number of postbuffer iteration to save
    dynModel                    = false     # boolean variable to tell the simulator which model to use (dynModel=True-->it'll use dynamic model, dynModel=False-->it'll use kinematic model)

    selectedStates.Np           = 30                            # Number of points to take from each previous trajectory to build the convex hull
    selectedStates.Nl           = 5                             # Number of previous laps to include in the convex hull
    Nl                          = selectedStates.Nl
    selectedStates.selStates    = zeros(Nl*selectedStates.Np,4)  
    selectedStates.statesCost   = zeros(Nl*selectedStates.Np)

    mpcParams.N                 = 10                        #lenght of prediction horizon
    mpcParams.vPathFollowing    = 0.6                       # reference velocity for the path following stage
    mpcParams.QderivZ           = 0.1*[0,0.0,0.1,0.1]       # weights for the states in the derivative cost
    mpcParams.QderivU           = 0.1*[1.0,10.0]            # weights for the control inputs in the derivative cost
    mpcParams.R                 = 0.0*[1.0,1.0]             # weights on the control inputs
    mpcParams.Q                 = [0.0,50.0,0.1,10.0]       # weights on the states for path following
    mpcParams.Q_cost            = 0.7                       # weight on the cost to get from a given point to the targe
    mpcParams.Q_lane            = 1.0                       # weight on the soft constraint on the Q_lane
    mpcParams.Q_alpha           = 1.0                       # weight on the soft constraint for convex hull
    mpcParams.Q_vel             = 1.0                       # weight on the soft constraint fot the max velocity
    mpcParams.Q_obs             = ones(Nl*selectedStates.Np)# weight to esclude some of the old trajectories
    mpcParams.Q_ell             = [0.001,0.001,0.024]#[0.01,0.01,0.24]          # weights defining the ellipse for the obstacle

    trackCoeff.nPolyCurvature   = 4                       # 4th order polynomial for curvature approximation
    trackCoeff.nPolyXY          = 6 
    trackCoeff.coeffCurvature   = [0.0;0.0;0.0;0.0;0.0]         # polynomial coefficients for curvature approximation (zeros for straight line)
    trackCoeff.width            = 2                   # width of the track (0.6m)
    trackCoeff.ds               = 1//10#4//100#1//10 # is defined as a rational number so we can use it to calculate indices in matrix. with float becomes error

    modelParams.l_A             = 0.125
    modelParams.l_B             = 0.125 
    modelParams.dt              = 0.1
    modelParams.u_lb            = ones(mpcParams.N,1) * [-0.6  -pi/3]  #-0.6 for braking    # lower bounds on control actions
    modelParams.u_ub            = ones(mpcParams.N,1) * [ 0.6   pi/3]       #1.2           # upper bounds on control actions
    modelParams.z_lb            = ones(mpcParams.N+1,1) * [-Inf -Inf -Inf  0]           # lower bounds on states
    modelParams.z_ub            = ones(mpcParams.N+1,1) * [ Inf  Inf  Inf  Inf]            # upper bounds on states
    modelParams.v_max           = 2.5                                                      # maxium velocity for trajectory tracking
    modelParams.max_alpha       = 10                                                       # maximum slip angle 
    modelParams.mass            = 1.98 # kg
    modelParams.mu              = 0.85
    modelParams.g               = 9.81 # m/s^2
    modelParams.I_z             = 0.03 # kg * m^2
    modelParams.B               = 6.0  
    modelParams.C               = 1.6  

    oldTraj.n_oldTraj           = simVariables.n_laps                                                     # number of old Trajectories for safe set
    oldTraj.oldTraj             = zeros(buffersize,4,oldTraj.n_oldTraj)                  # old trajectories in s-ey frame
    oldTraj.oldTrajXY           = zeros(buffersize,6,oldTraj.n_oldTraj)                  # old trajectories in x-y frame
    oldTraj.oldNIter            = zeros(oldTraj.n_oldTraj)                               # how many iterations to complete that particular iteration
    oldTraj.oldInput            = zeros(buffersize,2,oldTraj.n_oldTraj)                  # old Inputs 
    oldTraj.costs               = zeros(7,buffersize,oldTraj.n_oldTraj)                  # optimal values of the cost function elements for each optimization iteration
    oldTraj.z_pred_sol          = zeros(mpcParams.N+1,4,buffersize,oldTraj.n_oldTraj)    # predicted states (in s-ey frame) for each iteration of past rounds
    oldTraj.u_pred_sol          = zeros(mpcParams.N,2,buffersize,oldTraj.n_oldTraj)      # predicted input for each iteration of past rounds
    oldTraj.cost2target         = zeros(buffersize,oldTraj.n_oldTraj)                    # number of iterations needed to arrive at the target
    oldTraj.curvature           = zeros(buffersize,oldTraj.n_oldTraj)                    # all the curvatures calculated in each iteration of each lap
    oldTraj.oldAlpha            = zeros(Nl*selectedStates.Np,buffersize,oldTraj.n_oldTraj)# all the alphas from all iterations of all LMPC laps
    oldTraj.costLap             = zeros(oldTraj.n_oldTraj)                               # number of iterations to complete a full lap
    oldTraj.data_log            = zeros(trackCoeff.nPolyXY +1,3,buffersize,oldTraj.n_oldTraj)# logs all the data needed to perform the offline change of coordinates from s-ey to x-y



    mpcSol.u                    = zeros(mpcParams.N,2)              # array containing all the control inputs computed by the MPC at a given iteration
    mpcSol.z                    = zeros(mpcParams.N+1,4)            # array containing all the states computed by the MPC at a given iteration
    mpcSol.cost                 = zeros(5)                          # optimal costs as computed by the MPC at a given iteration
    mpcSol.alpha                = zeros(buffersize,Nl*selectedStates.Np,oldTraj.n_oldTraj)  # coefficients of the convex hull
    
    lapStatus.currentLap        = 1         # initialize lap number
    lapStatus.currentIt         = 0         # current iteration in lap 

    obstacle.obstacle_active    = false     # true if we have t consider the obstacles in the optimization problem
    obstacle.lap_active         = 9         # number of the first lap in which the obstacles are used
    obstacle.obs_detect         = 10         # maximum distance at which we can detect obstacles (in terms of s!!)
    obstacle.n_obs              = 1         # number of obstacles
    obstacle.s_obs_init         = [20]    # initial s coordinate of each obstacle
    obstacle.ey_obs_init        = [-0.8]       # initial ey coordinate of each obstacle
    obstacle.v_obs_init         = [0]       # initial velocity of each obstacles
    obstacle.r_s                = 0.5
    obstacle.r_ey               = 0.2
    obstacle.inv_step           = 3         # number of step of invariance required for the safe set

    
end