using JLD
using PyPlot
using PyCall
@pyimport matplotlib.animation as anim
@pyimport matplotlib.patches as patch
using JLD, ProfileView

matplotlib[:style][:use]("classic") # somehow my julia version changed plotting style 



include("../helper/classes.jl")
include("../helper/prova.jl")
include("createBorders.jl")
include("xyPredictions.jl")
include("xyObstacle.jl")


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
    obs_log    = Data["obs"]

    v_ref = mpcParams.vPathFollowing

    inner_x,inner_y,outer_x,outer_y = createBorders(x_track,y_track,trackCoeff,oldTraj)

    max_cost = findmax(oldTraj.costLap[:])[1]

    pred_sol_xy_1=xyObstacle(oldTraj,obs_log,1,laps[1],trackCoeff)
    pred_sol_xy_2=xyObstacle(oldTraj,obs_log,2,laps[1],trackCoeff)
    #   pred_sol_xy_3=xyObstacle(oldTraj,obs_log,3,laps[1],trackCoeff)

    close("all")

    figure()
    plot(oldTraj.costLap[:])
    grid("on")
    title("Cost of laps")

    vel_mean=zeros(laps[length(laps)])

    for j=1:laps[length(laps)]
        vel_mean[j] = mean(oldTraj.oldTraj[:,4,j])
    end

    figure()
    plot(vel_mean)
    plot(2.5*ones(laps[length(laps)]),"-")
    ylim(0,2.7)
    grid("on")
    title("Mean velocity over laps")

    for i = laps

        println("mean velocity of lap $i =",mean(oldTraj.oldTraj[:,4,i]))
        println("cost of lap $i= ",oldTraj.costLap[i])


        x = oldTraj.oldTrajXY[1:oldTraj.costLap[i],1,i]
        y = oldTraj.oldTrajXY[1:oldTraj.costLap[i],2,i]
        
        ellfig = figure()
        ax = ellfig[:add_subplot](1,1,1)
        ax[:set_aspect]("equal")
        plot(x,y,"g")
        plot(x_track',y_track',"r",inner_x,inner_y,"b",outer_x,outer_y,"b")
        ell1 = patch.Ellipse([pred_sol_xy_1[1,1],pred_sol_xy_1[2,1]], 1, 0.4, angle=0.0)
        ell2 = patch.Ellipse([pred_sol_xy_2[1,1],pred_sol_xy_2[2,1]], 1, 0.4, angle=90.0)
        #ell3 = patch.Ellipse([pred_sol_xy_3[1,1],pred_sol_xy_3[2,1]], 1, 0.4, angle=90.0)
        ax[:add_artist](ell1)
        ax[:add_artist](ell2)
        #ax[:add_artist](ell3)
        axis("equal")
        grid("on")
        title("X-Y view of Lap $i")

        #println("x= ",x)
        #println("y= ",y)

        t=linspace(1,max_cost,max_cost)

        figure()
        subplot(221)
        plot(t,oldTraj.oldTraj[1:max_cost,1,i])
        title("State S in lap $i ")
        grid("on")
        subplot(222)
        plot(t,oldTraj.oldTraj[1:max_cost,2,i]) 
        title("State ey in lap $i ")
        grid("on")
        subplot(223)
        plot(t,oldTraj.oldTraj[1:max_cost,3,i])
        title("State epsi in lap $i ")
        grid("on")
        subplot(224)
        plot(t,oldTraj.oldTraj[1:max_cost,4,i])
        legend(["State v"])
        title("State v in lap $i ")
        grid("on")

    end


end


function eval_pred(code::AbstractString,laps=Array{Int64})

    file = "../data/$(code).jld" 

    Data = load(file)

    x_track = Data["x_track"]
    y_track = Data["y_track"]
    trackCoeff = Data["trackCoeff"]
    modelParams = Data["modelParams"]
    mpcParams = Data["mpcParams"]
    buffersize = Data["buffersize"]
    oldTraj     = Data["oldTraj"]
    selectedStates = Data["selectedStates"]
    obs_log    = Data["obs"]

    v_ref = mpcParams.vPathFollowing

    inner_x,inner_y,outer_x,outer_y = createBorders(x_track,y_track,trackCoeff,oldTraj)




    t  = linspace(1,10,10)
    tN = linspace(1,11,11)

   

    for i = laps

        pred_sol_xy_obs=xyObstacle(oldTraj,obs_log,1,i,trackCoeff)
        pred_sol_xy = xyPredictions(oldTraj,i,trackCoeff)

        #println(obs_log[1:100,1,1,10])


        # s_pred    = oldTraj.z_pred_sol[:,1,1,i]
        # ey_pred   = oldTraj.z_pred_sol[:,2,1,i]
        # olds    = selectedStates[1:10,1,1,i]
        # olds2    = selectedStates[11:20,1,1,i]
        # oldey   = selectedStates[1:10,2,1,i]
        # oldey2   = selectedStates[11:20,2,1,i]
        
        # figure()
        # plot(s_pred,ey_pred,"or")
        # plot(olds,oldey,"b")
        # plot(olds2,oldey2,"b")
        # title("State S in lap $i ")
        # grid("on")


        for j = 2:2000

          
            #clf()
        
        # s    = oldTraj.oldTraj[:,1,i]
        # ey   = oldTraj.oldTraj[:,2,i]
        # epsi = oldTraj.oldTraj[:,3,i]
        # v    = oldTraj.oldTraj[:,4,i]

            s_pred    = oldTraj.z_pred_sol[:,1,j,i]
            ey_pred   = oldTraj.z_pred_sol[:,2,j,i]
            epsi_pred = oldTraj.z_pred_sol[:,3,j,i]
            v_pred    = oldTraj.z_pred_sol[:,4,j,i]

            olds     = selectedStates[1:20,1,j,i]
            olds2    = selectedStates[21:40,1,j,i]
            oldey    = selectedStates[1:20,2,j,i]
            oldey2   = selectedStates[21:40,2,j,i]
            oldepsi  = selectedStates[1:20,3,j,i]
            oldepsi2 = selectedStates[21:40,3,j,i]
            oldv     = selectedStates[1:20,4,j,i]
            oldv2    = selectedStates[21:40,4,j,i]

            x_pred  = pred_sol_xy[:,1,j]
            y_pred  = pred_sol_xy[:,2,j]


            
            
            figure(1)
            clf()
            subplot(221)
            plot(s_pred,ey_pred,"or")
            plot(olds,oldey,"b")
            plot(olds2,oldey2,"b")
            #ylim(findmin(oldTraj.z_pred_sol[:,2,:,i])[1],findmax(oldTraj.z_pred_sol[:,2,:,i])[1])
            title("State ey in lap $i ")
            grid("on")

            subplot(222)
            plot(s_pred,epsi_pred,"or")
            plot(olds,oldepsi,"b")
            plot(olds2,oldepsi2,"b")
            #ylim(findmin(oldTraj.z_pred_sol[:,3,:,i])[1],findmax(oldTraj.z_pred_sol[:,3,:,i])[1])
            title("State epsi in lap $i ")
            grid("on")

            subplot(223)
            plot(s_pred,v_pred,"or")
            plot(olds,oldv,"b")
            plot(olds2,oldv2,"b")
            #ylim(findmin(oldTraj.z_pred_sol[:,4,:,i])[1],findmax(oldTraj.z_pred_sol[:,4,:,i])[1])
            title("State v in lap $i ")
            grid("on")


            ellfig = figure(2)
            clf()
            ax = ellfig[:add_subplot](1,1,1)
            ax[:set_aspect]("equal")
            angle_ell = atan2(pred_sol_xy_obs[2,j]-(pred_sol_xy_obs[2,j-1]),pred_sol_xy_obs[1,j]-(pred_sol_xy_obs[1,j-1]))
            angle_deg = (angle_ell*180)/pi
            ell1 = patch.Ellipse([pred_sol_xy_obs[1,j],pred_sol_xy_obs[2,j]], 1, 0.4, angle=angle_deg)
            ax[:add_artist](ell1)
            x = oldTraj.oldTrajXY[j,1,i]
            y = oldTraj.oldTrajXY[j,2,i]
            x_obs = pred_sol_xy_obs[1,j]
            y_obs = pred_sol_xy_obs[2,j]
            plot(x_track',y_track',"g",inner_x,inner_y,"b",outer_x,outer_y,"b")
            plot(x,y,"or")
            plot(x_obs,y_obs,"ob")
            #println(pred_sol_xy_obs[1,j])
            

            sleep(0.0001)
        end
    end
end




function animation(code::AbstractString,lap=Int64)



    file = "../data/$(code).jld" 

    Data = load(file)

    x_track = Data["x_track"]
    y_track = Data["y_track"]
    trackCoeff = Data["trackCoeff"]
    modelParams = Data["modelParams"]
    mpcParams = Data["mpcParams"]
    buffersize = Data["buffersize"]
    oldTraj     = Data["oldTraj"]
    obs_log    = Data["obs"]

    N = mpcParams.N
    predTraj = oldTraj.z_pred_sol
    predInp  = oldTraj.u_pred_sol

    v_ref = mpcParams.vPathFollowing

    inner_x,inner_y,outer_x,outer_y = createBorders(x_track,y_track,trackCoeff,oldTraj)
    pred_sol_xy_1=xyObstacle(oldTraj,obs_log,1,lap,trackCoeff)
    


    #Construct Figure and Plot Data
    fig = figure("PathFollowing",figsize=(10,10))
    ax = axes(xlim = (-15,15),ylim=(-15,4))
    global line1 = ax[:plot]([],[],"or")[1]
    global line2 = ax[:plot]([],[],"b-")[1]
    global line3 = ax[:plot]([],[],"g-")[1]
    global line4 = ax[:plot]([],[],"b-")[1]
    ell1 = patch.Ellipse([pred_sol_xy_1[1,1],pred_sol_xy_1[2,1]], 1, 0.4, angle=0.0)
    
    global line5 = ax[:add_artist](ell1)
   
    



    function init()
        global line1
        global line2
        global line3
        global line4
        global line5
        
        line1[:set_data]([],[])
        line2[:set_data]([],[])
        line3[:set_data]([],[])
        line4[:set_data]([],[])
    
        return (line1,line2,line3,line4,line5,Union{})  # Union{} is the new word for None
    end

    steps = 15

    function animate(i)
        k = i+2
        global line1
        global line2
        global line3
        global line4
        global line5
        
        
        x = oldTraj.oldTrajXY[k,1,lap]
        y = oldTraj.oldTrajXY[k,2,lap]
        angle_ell = atan2(pred_sol_xy_1[2,k]-(pred_sol_xy_1[2,k-1]),pred_sol_xy_1[1,k]-(pred_sol_xy_1[1,k-1]))
        angle_deg = (angle_ell*180)/pi
        ell1 = patch.Ellipse([pred_sol_xy_1[1,k],pred_sol_xy_1[2,k]], 1, 0.4, angle=angle_deg)
        
        
        line1[:set_data](x,y)
        line2[:set_data](inner_x[:],inner_y[:])
        line3[:set_data](x_track[:]',y_track[:]')
        line4[:set_data](outer_x[:],outer_y[:])     
        line5 = ax[:add_artist](ell1)
           
        return (line1,line2,line3,line4,line5,Union{})
    end

    myanim = anim.FuncAnimation(fig, animate, init_func=init, frames=1000, interval=20)

    myanim[:save]("3Lines.mp4", bitrate=-1, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])


    function html_video(filename)
        base64_video = base64encode(open(readbytes, filename))
        """<video controls src="data:video/x-m4v;base64,$base64_video">"""
    end

end

