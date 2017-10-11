using JLD
using PyPlot
using PyCall
@pyimport matplotlib.animation as anim
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
        plot(x,y,"g")
        plot(x_track',y_track',"r",inner_x,inner_y,"b",outer_x,outer_y,"b")
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


function animation(code::AbstractString,laps=Array{Int64})



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

    #Construct Figure and Plot Data
    fig = figure("PathFollowing",figsize=(10,10))
    ax = axes(xlim = (-15,15),ylim=(-15,4))
    global line1 = ax[:plot]([],[],"or")[1]
    global line2 = ax[:plot]([],[],"b-")[1]
    global line3 = ax[:plot]([],[],"g-")[1]
    global line4 = ax[:plot]([],[],"b-")[1]


    function init()
        global line1
        global line2
        global line3
        global line4
        line1[:set_data]([],[])
        line2[:set_data]([],[])
        line3[:set_data]([],[])
        line4[:set_data]([],[])
    
        return (line1,line2,line3,line4,Union{})  # Union{} is the new word for None
    end

    steps = 15

    function animate(i)
        k = i+1
        global line1
        global line2
        global line3
        global line4
        
        x = oldTraj.oldTrajXY[k,1,1]
        y = oldTraj.oldTrajXY[k,2,1]
        
        line1[:set_data](x,y)
        line2[:set_data](inner_x[:],inner_y[:])
        line3[:set_data](x_track[:]',y_track[:]')
        line4[:set_data](outer_x[:],outer_y[:])        
        return (line1,line2,line3,line4,Union{})
    end

    myanim = anim.FuncAnimation(fig, animate, init_func=init, frames=1000, interval=20)

    myanim[:save]("3Lines.mp4", bitrate=-1, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])


    function html_video(filename)
        base64_video = base64encode(open(readbytes, filename))
        """<video controls src="data:video/x-m4v;base64,$base64_video">"""
    end

end

