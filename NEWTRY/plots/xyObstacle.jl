# This file translates the obstacles position from s-ey to x-y frame (for plotting purposes)


function xyObstacle(oldTraj::classes.OldTrajectory,obs_log::Array{Float64},obstacle::Int64,lap::Int64,trackCoeff::classes.TrackCoeff)

    obs   = obs_log[:,:,obstacle,lap]

    

    buffersize = size(obs)[1]
    

    data_log   = oldTraj.data_log

    OrderXY        = trackCoeff.nPolyXY
    OrderThetaCurv = trackCoeff.nPolyCurvature

    n_poly = 41

    ds = 1//10

    s_vec = zeros(OrderXY+1)

    pred_sol_xy = zeros(2,buffersize,1)

    include("../helper/loadTestMap.jl")
    x_track, y_track = loadTestMap() #load a track to test controller

    println("x_track= ",size(x_track))
    println("nodes= ",size([x_track; y_track]))
    for i = 1:buffersize

            nodes          = [x_track; y_track]
            n_nodes        = size(x_track)[2]


            s_nearest = obs[i,1]

            idx_min = floor(s_nearest)*10


            idx_start = idx_min - 10
            idx_end   = idx_min + 30

            s_start = (idx_start-1)*ds



            if idx_start<=0
                 nodes_XY = hcat(nodes[:,n_nodes+idx_start:n_nodes],nodes[:,1:idx_end])       # then stack the end and beginning of a lap together
            #     #nodes_Y = hcat(nodes[2,n_nodes+idx_start:n_nodes],nodes[2,1:idx_end])
               
                 #idx_start = n_nodes+idx_start
            elseif idx_end>=n_nodes                   # if the end is behind the finish line
                 nodes_XY = hcat(nodes[:,idx_start:n_nodes],nodes[:,1:idx_end-n_nodes])       # then stack the end and beginning of the lap together
                 #nodes_Y = hcat(nodes[2,idx_start:n_nodes],nodes[2,1:idx_end-n_nodes])
            else                               # if we are somewhere in the middle of the track
                nodes_XY = nodes[:,idx_start:idx_end]     # then just use the nodes from idx_start to end for interpolation
             #nodes_Y = nodes[2,idx_start:idx_end]
            end


            nodes_X = vec(nodes_XY[1,:])
            nodes_Y = vec(nodes_XY[2,:])


            itp_matrix = zeros(n_poly,OrderXY+1)

            for ind=1:n_poly
                for k=0:OrderXY

                    itp_matrix[ind,OrderXY+1-k] = (s_start + (ind-1)*ds)^k
                end
            end



            coeffY = itp_matrix\nodes_Y
            coeffX = itp_matrix\nodes_X

            
          
            s_vec = zeros(OrderXY+1)::Array{Float64}
            sdot_vec = zeros(OrderXY+1)::Array{Float64}

            for k = 1:OrderXY+1
                    s_vec[k] = obs[i,1]^(OrderXY-k+1)
                    
            end
            for k = 1:OrderXY
                    sdot_vec[k] = (OrderXY+1-k)* obs[i,1]^(OrderXY-k)
            end


            XCurve  = dot(coeffX,s_vec)
            YCurve  = dot(coeffY,s_vec)

            dX = dot(coeffX,sdot_vec)
            dY = dot(coeffY,sdot_vec)      

            
            xyPathAngle = atan2(dY,dX)

            pred_sol_xy[2,i] = YCurve + obs[i,2]*cos(xyPathAngle)
            pred_sol_xy[1,i] = XCurve - obs[i,2]*sin(xyPathAngle)

        
    end

    return pred_sol_xy
end

