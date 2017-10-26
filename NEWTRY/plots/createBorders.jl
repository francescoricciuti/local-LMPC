# This file creates the track borders

using JuMP
using Ipopt
using JLD
function createBorders(x_track::Array{Float64},y_track::Array{Float64},trackCoeff::classes.TrackCoeff,oldTraj::classes.OldTrajectory)


    states=zeros(1,6)
    inner_x = zeros(size(x_track)[2])
    inner_y = zeros(size(x_track)[2])
    outer_x = zeros(size(x_track)[2])
    outer_y = zeros(size(x_track)[2])

    width   = trackCoeff.width

    for i=1:size(x_track)[2]
       
        x_now = x_track[i]
        y_now = y_track[i]
        states[1,1] = x_now
        states[1,2] = y_now

        zCurr_s, coeffCurv, XCurve, YCurve, xyPathAngle = trackFrameConversion(states,x_track,y_track,trackCoeff, 1,oldTraj,1,false)

        inner_y[i] = YCurve + (-width/2)*cos(xyPathAngle)
        inner_x[i] = XCurve - (-width/2)*sin(xyPathAngle)
        outer_y[i] = YCurve + (width/2)*cos(xyPathAngle)
        outer_x[i] = XCurve - (width/2)*sin(xyPathAngle)
    end


    return inner_x,inner_y,outer_x,outer_y
end