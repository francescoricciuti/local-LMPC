# This file defines all the Julia Types needed by the problem

module classes# VARIOUS TYPES FOR CALCULATIONS
using JuMP

type LapStatus
    currentLap::Int64  # current lap number
    currentIt::Int64   # current iteration in the current lap

    LapStatus(currentLap = 1, currentIt = 1) = new(currentLap,currentIt)
end

type SimulationVariables
    buffersize::Int64  # used to initialize the dimensions of the variables in which we will save the data of the Simulations 
    n_laps::Int64      # number of laps we want to simulate 
    n_pf::Int64        # number of path following laps (must be at least 2)
    postbuff::Int64    # number of postbuffer iteration to save
    dynModel::Bool     # boolean variable to tell the simulator which model to use (dynModel=True-->it'll use dynamic model, dynModel=False-->it'll use kinematic model)

    SimulationVariables(buffersize = 2000, n_laps = 30, n_pf = 3,postbuff=30,dynModel=true) = new(buffersize,n_laps,n_pf,postbuff,dynModel)

end


type OldTrajectory                  # informations about previous trajectories
    n_oldTraj::Int64                # number of Old Trajectories
    oldTraj::Array{Float64}         # Old Trajectories in coordinates relative to the race track
    oldTrajXY::Array{Float64}       # Old Trajectories in XY coordinates
    oldNIter::Array{Float64}        # number of iteration needed to complete that particular lap
    oldInput::Array{Float64}        # Old inputs
    costs::Array{Float64}           # Old values of cost functions of old optimization rounds
    z_pred_sol::Array{Float64}      # all the states predicted by the MPC
    u_pred_sol::Array{Float64}      # all the control actions predicted by the MPC
    cost2target::Array{Float64}     # cost to arrive at the target, i.e. how many iterations from the start to the end of the lap
    curvature::Array{Float64}       # all the curvature calculated in each step of each lap
    oldAlpha::Array{Float64}        # all the alphas computed in each iteration of each LMPC lap
    costLap::Array{Int64}           # number of iterations to complete a full lap
    data_log::Array{Float64}        # logs all the data needed to perform the offline change of coordinates from s-ey to x-y
    
   
    OldTrajectory(n_oldTraj = 0, oldTraj=Float64[],oldTrajXY=Float64[],oldNIter=Float64[],oldInput=Float64[],costs=Float64[],z_pred_sol=Float64[],
                  u_pred_sol=Float64[],cost2target= Float64[],curvature=Float64[],oldAlpha=Float64[],costLap=Int64[],data_log=Float64[]) =
                 new(n_oldTraj, oldTraj,oldTrajXY,oldNIter,oldInput,costs,z_pred_sol,u_pred_sol,cost2target,curvature,oldAlpha,costLap,data_log)
end

type SelectedStates                 # Values needed for the convex hull formulation
    selStates::Array{Float64}       # selected states from previous laps ...
    statesCost::Array{Float64}      # ... and their related costs
    Np::Int64                       # number of states to select from each previous lap
    SelectedStates(selStates=Float64[],statesCost=Float64[],Np=6) = new(selStates,statesCost,Np)
end


type MpcParams                     # parameters for MPC 
    N::Int64                       # prediction horizon
    vPathFollowing::Float64        # reference velocity for the path following stage
    QderivZ::Array{Float64,1}      # weights on the states for the derivative cost
    QderivU::Array{Float64,1}      # weights on the control inputs for the derivative cost
    R::Array{Float64,1}            # weights on the control inputs
    Q::Array{Float64,1}            # weights on the states for path following
    Q_cost::Float64                # weight on the cost to get from a given point to the target
    Q_lane::Float64                # weight on the soft constraint for the lane
    Q_alpha::Float64               # weight on the soft constraint for the convex hull

    MpcParams(N=0,vPathFollowing=1.0,QderivZ=Float64[],QderivU=Float64[],R=Float64[],Q=Float64[],Q_cost=0.7,Q_lane=0.5,Q_alpha=1.0) = 
    new(N,vPathFollowing,QderivZ,QderivU,R,Q,Q_cost,Q_lane,Q_alpha)
end

type PosInfo            # current position information
    s_start::Float64    # first point with s>0
    s::Float64          # current point
    s_target::Float64   # s coordinate of the target
    PosInfo(s_start=0.0,s=0.0,s_target=0.0) = new(s_start,s,s_target)
end

type MpcSol                      # MPC solution output
    a_x::Float64                 # accelleration computed by the MPC (only the first one, i.e. the one that will be applied to the system)
    d_f::Float64                 # steering computed by the MPC (only the first one, i.e. the one that will be applied to the system)
    solverStatus::Symbol         # status of the solver
    u::Array{Float64}            # array containing all the inputs computed by the MPC
    z::Array{Float64}            # array containing all the inputs computed by the MPC
    cost::Array{Float64}         # array containing the values of the cost functions of this particular MPC optimization
    alpha::Array{Float64}        # parameters of the convex hull
    MpcSol(a_x=0.0,d_f=0.0,solverStatus=Symbol(),u=Float64[],z=Float64[],cost=Float64[],alpha=Float64[]) = new(a_x,d_f,solverStatus,u,z,cost,alpha)
end


type TrackCoeff         # coefficients of track
    coeffAngle::Array{Float64,1}
    coeffCurvature::Array{Float64,1}
    nPolyCurvature::Int64               # order of the interpolation polynom
    nPolyXY::Int64                      # order of the interpolation polynom of the x y coordinates
    width::Float64                      # lane width -> is used in cost function as soft constraints (to stay on track)
    ds::Float64
    TrackCoeff(coeffAngle=Float64[], coeffCurvature=Float64[], nPolyCurvature=4, nPolyXY = 6, width=1.0, ds=0.1) = 
    new(coeffAngle,coeffCurvature,nPolyCurvature,nPolyXY,width,ds)
end

type ModelParams                # Values for the model
    l_A::Float64                # distance from CoM to front wheel
    l_B::Float64                # distance from CoM to rear wheel
    dt::Float64                 # time step
    u_lb::Array{Float64}        # lower bounds for control inputs
    u_ub::Array{Float64}        # upper bounds for control inputs
    z_lb::Array{Float64}        # lower bounds for states
    z_ub::Array{Float64}        # upper bounds for states
    v_max::Float64              # maximum velocity for path following
    max_alpha::Float64          # maximum slip angle allowed
    mass::Float64              
    mu::Float64
    g::Float64
    I_z::Float64
    B::Float64
    C::Float64
    ModelParams(l_A=0.25,l_B=0.25,dt=0.1,u_lb=Float64[],u_ub=Float64[],z_lb=Float64[],z_ub=Float64[],v_max = 2.0,max_alpha = 10, mass = 1.98, mu = 0.85, g = 9.81, I_z = 0.03, B =3.0, C = 1.25) = 
    new(l_A,l_B,dt,u_lb,u_ub,z_lb,z_ub,v_max,max_alpha, mass, mu, g, I_z, B, C)
end

end

