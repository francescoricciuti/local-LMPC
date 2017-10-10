using JLD
using PyPlot
using PyCall
@pyimport matplotlib.animation as animation
using JLD, ProfileView

matplotlib[:style][:use]("classic") # somehow my julia version changed plotting style 



include("../helper/classes.jl")
include("../helper/prova.jl")
include("createBorders.jl")


function eval_sim(code::AbstractString,laps=Array{Int64})

    file = "../data/$(code).jld" 

    Data = load(file)

    x_track = Data["x_track"]
    y_track = Data["y_track"]
    trackCoeff = Data["trackCoeff"]
    modelParams = Data["modelParams"]
    mpcParams = Data["mpcParams"]
    buffersize = Data["buffersize"]
    oldTraj     = Data["oldTraj"]

    v_ref = mpcParams.vPathFollowing

    inner_x,inner_y,outer_x,outer_y = createBorders(x_track,y_track,trackCoeff)



    close("all")

    for i = laps

        println("mean velocity of lap $i =",mean(oldTraj.oldTraj[:,4,i]))

        x = oldTraj.oldTrajXY[:,1,i]
        y = oldTraj.oldTrajXY[:,2,i]
    
        figure()
        plot(x,y)
        plot(x_track',y_track',inner_x,inner_y,outer_x,outer_y)
        axis("equal")
        grid("on")
        title("X-Y view of Lap $i")

        t=linspace(1,size(oldTraj.oldTraj[:,2,i])[1],size(oldTraj.oldTraj[:,2,i])[1])

        figure()
        subplot(221)
        plot(t,oldTraj.oldTraj[:,1,i])
        title("State S in lap $i ")
        grid("on")
        subplot(222)
        plot(t,oldTraj.oldTraj[:,2,i])
        title("State ey in lap $i ")
        grid("on")
        subplot(223)
        plot(t,oldTraj.oldTraj[:,3,i])
        title("State epsi in lap $i ")
        grid("on")
        subplot(224)
        plot(t,oldTraj.oldTraj[:,4,i])
        legend(["State v"])
        title("State v in lap $i ")
        grid("on")

    end


end

function eval_pred(code::AbstractString,laps=Array{Int64},iter=Array{Int64})

    file = "../data/$(code).jld" 

    Data = load(file)

    x_track     = Data["x_track"]
    y_track     = Data["y_track"]
    trackCoeff  = Data["trackCoeff"]
    modelParams = Data["modelParams"]
    mpcParams   = Data["mpcParams"]
    buffersize  = Data["buffersize"]
    oldTraj     = Data["oldTraj"]

    close("all")
    println("starting states= ",oldTraj.oldTraj[2,:,1])
    println("predicted states= ",oldTraj.z_pred_sol[2,:,2,1])
    println("next states= ",oldTraj.oldTraj[3,:,1])
    println("predicted  inputs= ",oldTraj.u_pred_sol[2,:,2,1])

    # for i = laps
    #     for k = iter
    #         s_pred = oldTraj.z_pred_sol[:,1,k,i]

    #         println(oldTraj.z_pred_sol[:,2,k,i])

    #         ey_pred = oldTraj.z_pred_sol[:,2,k,i]

    #         epsi_pred = oldTraj.z_pred_sol[:,3,k,i]

    #         v_pred = oldTraj.z_pred_sol[:,4,k,i]

    #         t=linspace(1,mpcParams.N+1,mpcParams.N+1)

    #         figure()
    #         plot(t,s_pred)
    #         title("s, Lap $i, iteration $k ")
    #         grid("on")


    #         figure()
    #         plot(t,ey_pred)
    #         title("ey, Lap $i, iteration $k ")
    #         grid("on")


    #         figure()
    #         plot(t,epsi_pred)
    #         title("epsi, Lap $i, iteration $k ")
    #         grid("on")


    #         figure()
    #         plot(t,v_pred)
    #         title("v, Lap $i, iteration $k ")
    #         grid("on")

            
    
    #     end
    # end


end   