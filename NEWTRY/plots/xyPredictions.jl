# This file translates the states predicted by the MPC from s-ey to x-y frame (for plotting purposes)

function xyPredictions(oldTraj::classes.OldTrajectory,lap::Int64,trackCoeff::classes.TrackCoeff)

    pred_sol   = oldTraj.z_pred_sol[:,:,:,lap]

    buffersize = size(pred_sol)[3]
    N          = size(pred_sol)[1]

    data_log   = oldTraj.data_log

    OrderXY        = trackCoeff.nPolyXY
    OrderThetaCurv = trackCoeff.nPolyCurvature

    n_poly = 41

    ds = 1//10

    s_vec = zeros(OrderXY+1)

    pred_sol_xy = zeros(N,2,buffersize,1)

    include("../helper/loadTestMap.jl")
    x_track, y_track = loadTestMap() #load a track to test controller

    
    for i = 1:buffersize


        for j = 1:N

            nodes          = [x_track; y_track]
            n_nodes        = size(x_track)[2]

            s_start   = (pred_sol[j,1,i] - 1)
            s_end     = (pred_sol[j,1,i] + 3)
            s_nearest = pred_sol[j,1,i]

            idx_start = 10*(floor(pred_sol[j,1,i]) - 1) 
            idx_end   = 10*(floor(pred_sol[j,1,i]) + 3)

            if idx_start>n_nodes
              idx_start=idx_start%n_nodes
              idx_end=idx_end%n_nodes
            end

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

            itp_matrix_curv = zeros(n_poly,OrderThetaCurv+1)

            for ind=1:n_poly
                for k=0:OrderThetaCurv

                    itp_matrix_curv[ind,OrderThetaCurv+1-k] = (s_start + (ind-1)*ds)^k
                end
            end
            
           # println("size of nodes x= ",size(nodes_X))
           # println("size of itpmatrix= ",size(itp_matrix))
           # println("s start= ",s_start)
           # println("s end= ",s_end)

            coeffY = itp_matrix\nodes_Y
            coeffX = itp_matrix\nodes_X
           

            b_curvature_vector = zeros(n_poly)

            Counter = 1
            

            for ind = 0:n_poly-1
                s_expression_der  = zeros(OrderXY+1)
                s_expression_2der = zeros(OrderXY+1)
                s_poly       = s_start + ind*ds
                for k=0:OrderXY-1
                    s_expression_der[OrderXY-k] = (k+1)*s_poly^k
                end
                for k=0:OrderXY-2
                    s_expression_2der[OrderXY-1-k] = (2+k*(3+k))*s_poly^k
                end

                dX  = dot(coeffX,s_expression_der)
                dY  = dot(coeffY,s_expression_der)
                ddX = dot(coeffX,s_expression_2der)
                ddY = dot(coeffY,s_expression_2der)

                curvature = (dX*ddY-dY*ddX)/((dX^2+dY^2)^(3/2)) #standard curvature formula

                b_curvature_vector[Counter] = curvature

                Counter = Counter + 1
            end


            
            coeffCurv  = itp_matrix_curv\b_curvature_vector

            s0 =  pred_sol[j,1,i]+0.001
          
            s_vec = zeros(OrderXY+1)::Array{Float64}
            sdot_vec = zeros(OrderXY+1)::Array{Float64}

            for k = 1:OrderXY+1
                    s_vec[k] = pred_sol[j,1,i]^(OrderXY-k+1)
                    
            end
            for k = 1:OrderXY
                    sdot_vec[k] = (OrderXY+1-k)* pred_sol[j,1,i]^(OrderXY-k)
            end


            XCurve  = dot(coeffX,s_vec)
            YCurve  = dot(coeffY,s_vec)

            dX = dot(coeffX,sdot_vec)
            dY = dot(coeffY,sdot_vec)      

            
            xyPathAngle = atan2(dY,dX)

            pred_sol_xy[j,2,i] = YCurve + pred_sol[j,2,i]*cos(xyPathAngle)
            pred_sol_xy[j,1,i] = XCurve - pred_sol[j,2,i]*sin(xyPathAngle)


            # current_s  = pred_sol[j,1,i]
            # current_ey = pred_sol[j,2,i]

            # for k=0:OrderXY
            #     s_vec[OrderXY-k+1]=current_s^k
            # end

            # coeffX = data_log[:,1,i,lap]
            # coeffY = data_log[:,2,i,lap]

           

            # pred_sol_xy[j,1,i] = dot(coeffX,s_vec) + current_ey*cos(data_log[1,3,i,lap])
            # pred_sol_xy[j,2,i] = dot(coeffY,s_vec) - current_ey*sin(data_log[1,3,i,lap])
        end
    end

    return pred_sol_xy
end

