module classes# VARIOUS TYPES FOR CALCULATIONS
using JuMP



type LapStatus
    currentLap::Int64       # current lap number
    currentIt::Int64        # current iteration in current lap
end

type MpcCoeff           # coefficients for trajectory approximation

    order::Int64
    pLength::Int64      # small values here may lead to numerical problems since the functions are only approximated in a short horizon
                        # "small" values are about 2*N, good values about 4*N
                        # numerical problems occur at the edges (s=0, when v is almost 0 and s does not change fast and at s=s_target)
    MpcCoeff( order=4, pLength=0) = new(order, pLength)
end



type OldTrajectory      # information about previous trajectories
    n_oldTraj::Int64                # number of Old Trajectories
    oldTraj::Array{Float64}         # Old Trajectories in coordinates relative to the race track
    oldTrajXY::Array{Float64}       # Old Trajectories in XY coordinates
    curvature::Array{Float64}       # curvature of the race track
    oldInput::Array{Float64}        # Old inputs
    oldNIter
    costs::Array{Float64}           # Old costs
    #lambda_sol::Array{Float64}
    z_pred_sol::Array{Float64}
    u_pred_sol::Array{Float64}
    #ssInfOn_sol::Array{Float64}
    #eps::Array{Float64}
    cost2Target::Array{Float64}    
    #copyInfo::Array{Float64}
    OldTrajectory(n_oldTraj = 0, oldTraj=Float64[],oldTrajXY=Float64[],curvature=Float64[],oldInput=Float64[],oldNIter=Float64[],
        costs=Float64[],z_pred_sol=Float64[],u_pred_sol=Float64[],cost2Target= Float64[]) =
                 new(n_oldTraj, oldTraj,oldTrajXY,curvature,oldInput,oldNIter,costs,z_pred_sol,u_pred_sol,cost2Target)
end

type SelectedStates
    selStates::Array{Float64}  # selected states from previous laps ...
    statesCost::Array{Float64} # ... and their related costs
    Np::Int64 # number of states to select from each previous lap
    SelectedStates(selStates=Float64[],statesCost=Float64[],Np=6) = new(selStates,statesCost,Np)
end



type MpcParams          # parameters for MPC solver
    N::Int64
    nz::Int64
    OrderCostCons::Int64
    Q::Array{Float64,1}
    Q_term::Array{Float64,1}
    Q_cost::Float64
    Q_lane::Float64
    Q_velocity::Float64
    R::Array{Float64,1}
    vPathFollowing::Float64
    QderivZ::Array{Float64,1}
    QderivU::Array{Float64,1}
    MpcParams(N=0,nz=0,OrderCostCons=0,Q=Float64[],Q_term=Float64[],Q_cost=1.0,Q_lane = 1.0,Q_velocity=1.0, R=Float64[],vPathFollowing=1.0,QderivZ=Float64[],QderivU=Float64[]) = new(N,nz,OrderCostCons,Q,Q_term,Q_cost,Q_lane, Q_velocity, R,vPathFollowing)
end

type PosInfo            # current position information
    s_start::Float64    # first point with s>0
    s::Float64          # current point
    s_target::Float64   # s coordinate of the target
    PosInfo(s_start=0.0,s=0.0,s_target=0.0) = new(s_start,s,s_target)
end

type MpcSol                      # MPC solution output
    a_x::Float64                 # accelleration computed by the MPC
    d_f::Float64                 # steering computed by the MPC
    solverStatus::Symbol         # status of the solver
    u::Array{Float64}            # array containing all the inputs computed by the MPC
    z::Array{Float64}            # array containing all the inputs computed by the MPC
    #lambda::Array{Float64,1}
    #ssInfOn::Array{Int64,1}
    #eps::Array{Float64}
    cost::Array{Float64}         # array containing the values of the cost functions
    MpcSol(a_x=0.0,d_f=0.0,solverStatus=Symbol(),u=Float64[],z=Float64[],cost=Float64[]) = new(a_x,d_f,solverStatus,u,z,cost)
end


type TrackCoeff         # coefficients of track
    coeffAngle::Array{Float64,1}
    coeffCurvature::Array{Float64,1}
    nPolyCurvature::Int64               # order of the interpolation polynom
    nPolyXY::Int64                      # order of the interpolation polynom of the x y coordinates
    width::Float64                      # lane width -> is used in cost function as soft constraints (to stay on track)
    ds::Rational{Int}
    TrackCoeff(coeffAngle=Float64[], coeffCurvature=Float64[], nPolyCurvature=4, nPolyXY = 6, width=1.0, ds=1/10) = new(coeffAngle,coeffCurvature,nPolyCurvature,nPolyXY,ds)
end

type ModelParams
    l_A::Float64                # distance from CoM to front wheel
    l_B::Float64                # distance from CoM to rear wheel
    dt::Float64                 # time step
    u_lb::Array{Float64}        # lower bounds for control inputs
    u_ub::Array{Float64}        # upper bounds for control inputs
    z_lb::Array{Float64}        # lower bounds for states
    z_ub::Array{Float64}        # upper bounds for states
    v_max::Float64              # maximum velocity for path following
    max_alpha::Float64
    mass::Float64
    mu::Float64
    g::Float64
    I_z::Float64
    B::Float64
    C::Float64
    ModelParams(l_A=0.25,l_B=0.25,dt=0.1,u_lb=Float64[],u_ub=Float64[],z_lb=Float64[],z_ub=Float64[],v_max = 2.0, max_alpha = 10.0, mass = 1.98, mu = 0.85, g = 9.81, I_z = 0.03, B =3.0, C = 1.25) = 
    new(l_A,l_B,dt,u_lb,u_ub,z_lb,z_ub,v_max, max_alpha, mass, mu, g, I_z, B, C)
end


end
