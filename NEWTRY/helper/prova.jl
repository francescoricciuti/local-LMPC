# This file performs the change of coordinates between x-y frame and s_ey frame

function trackFrameConversion(states_x::Array{Float64},x_track::Array{Float64},y_track::Array{Float64},trackCoeff::classes.TrackCoeff, itercount::Int64)

    OrderXY        = trackCoeff.nPolyXY
    OrderThetaCurv = trackCoeff.nPolyCurvature
    ds             = trackCoeff.ds

    x              = states_x[1]
    y              = states_x[2]
    psi            = states_x[5]


    nodes          = [x_track; y_track]
    n_nodes        = size(x_track)[2]


    #println("x_track= ",x_track[1:3])
    
    N_nodes_poly_back  = 10                                         
    N_nodes_poly_front = 30
    n_poly             = N_nodes_poly_front + N_nodes_poly_back +1  # overall number of points for interpolation

    dist_vec           = sqrt((x-x_track).^2 + (y-y_track).^2)
    dist_lane, idx_min = findmin(dist_vec)


    idx_start = idx_min - N_nodes_poly_back
    idx_end   = idx_min + N_nodes_poly_front

    #println("nodes= ",nodes)
    
    
    if idx_start<=0
        nodes_XY = hcat(nodes[:,n_nodes+idx_start:n_nodes],nodes[:,1:idx_end])       # then stack the end and beginning of a lap together
        #nodes_Y = hcat(nodes[2,n_nodes+idx_start:n_nodes],nodes[2,1:idx_end])
       
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
    #println("size nodes_x= ",size(nodes_X))
    #println("n_poly= ",n_poly)

    s_start   = (idx_start-1)*ds
    s_nearest = (idx_min-1)*ds
    s_nearest = round(s_nearest,4)

    
    #println("s start= ", s_start)
    #println("s nearest= ",s_nearest)


    # Create Matrix for interpolation
        # x-y-Matrix
        # The x-y-matrix is just filled with s^n values (n from 0 to the polynomial degree) for s values between s = 0 and s = n_poly*ds
        # with n_poly = number of points that approximate the polynomial

    itp_matrix = zeros(n_poly,OrderXY+1)

    for i=1:n_poly
        for k=0:OrderXY

            itp_matrix[i,OrderXY+1-k] = (s_start + (i-1)*ds)^k
        end
    end

    #println("itp_matrix= ",itp_matrix)


    itp_matrix_curv = zeros(n_poly,OrderThetaCurv+1)

    
    for i=1:n_poly
        for k=0:OrderThetaCurv

            itp_matrix_curv[i,OrderThetaCurv+1-k] = (s_start + (i-1)*ds)^k
        end
    end
    
   

    coeffY = itp_matrix\nodes_Y
    coeffX = itp_matrix\nodes_X
    # println("nodes_x= ",nodes_X)
    # println("coeffX= ",coeffX)
    # println("test= ",itp_matrix[11,:]*coeffX)

    b_curvature_vector = zeros(n_poly)

    Counter = 1
    

    for j = 0:n_poly-1
        s_expression_der  = zeros(OrderXY+1)
        s_expression_2der = zeros(OrderXY+1)
        s_poly       = s_start + j*ds
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

    #################################################

    discretization = 0.001 
    Counter = 1
    j_vec = zeros(OrderXY+1) 

    #println("s_nearest-ds= ",s_nearest-ds)
# if s_nearest >= ds
    dummy = size((s_nearest-ds):discretization:(s_nearest+ds))
    S_Value = zeros(dummy)  
    DistanceNew = zeros(dummy)
    # S_Value = zeros(0:discretization:2*ds)  
    # DistanceNew = zeros(0:discretization:2*ds)
    for j=(s_nearest-ds):discretization:(s_nearest+ds) 
        
        for i = 1:OrderXY+1
                j_vec[i] =j^(OrderXY-i+1)
        end


        XCurve = dot(coeffX, j_vec)
        YCurve = dot(coeffY, j_vec)

            #T just test
            # XCurve_t[Counter] = dot(coeffX, j_vec)
            # YCurve_t[Counter] = dot(coeffY, j_vec)

        S_Value[Counter] = j
        DistanceNew[Counter] = sqrt((x-XCurve).^2+(y-YCurve).^2) #distance of vehicle to every interpolated node
        Counter = Counter + 1
    end
    # println("size of j= ",size((s_nearest-ds):discretization:(s_nearest+ds)))
    # println("s_nearest= ",s_nearest)
    # println("ds= ",ds)
    # println("discretization= ",discretization)

    
    # else # just active while the car is nearest to the starting point as we cannot evaluate functon for negative s
    #     S_Value = zeros(0:discretization:1*ds) 
    #     DistanceNew = zeros(0:discretization:1*ds)
    #     for j=0:discretization:(ds) 
            
    #         for i = 1:OrderXY+1
    #             j_vec[i] =j^(OrderXY-i+1)
    #         end


    #         XCurve = dot(coeffX, j_vec)
    #         YCurve = dot(coeffY, j_vec)
    #         #println(j)
    #         #println(XCurve)
    #         #println(YCurve)

    #         #T just test
    #         # XCurve_t[Counter] = dot(coeffX, j_vec)
    #         # YCurve_t[Counter] = dot(coeffY, j_vec)

    #         S_Value[Counter] = j
    #         DistanceNew[Counter] = sqrt((x-XCurve).^2+(y-YCurve).^2) #distance of vehicle to every interpolated node
    #         Counter = Counter + 1
    #     end
    # end



    eyabs, idx_min_Dist = findmin(DistanceNew)

    println("eyabs= ",eyabs)

    s = S_Value[idx_min_Dist]

    # if s >= discretization

    #     s0 = s - discretization

    #     s0_vec = zeros(OrderXY+1)::Array{Float64}
    #     s_vec = zeros(OrderXY+1)::Array{Float64}
    #     sdot_vec = zeros(OrderXY+1)::Array{Float64}

    #     for i = 1:OrderXY+1
    #             s_vec[i] =s^(OrderXY-i+1)
    #             s0_vec[i] =s0^(OrderXY-i+1)
    #     end
    #     for i = 1:OrderXY
    #             sdot_vec[i] = (OrderXY+1-i)*s^(OrderXY-i)
    #     end

    #     XCurve0 = dot(coeffX,s0_vec)
    #     YCurve0 = dot(coeffY,s0_vec)

    #     XCurve  = dot(coeffX,s_vec)
    #     YCurve  = dot(coeffY,s_vec)

    #     dX = dot(coeffX,sdot_vec)
    #     dY = dot(coeffY,sdot_vec)      

    #     xyVectorAngle = atan2(y-YCurve0,x-XCurve0)
    #     xyPathAngle = atan2(dY,dX)

    # else
       
    s0 = s+discretization
    s0_vec = zeros(OrderXY+1)::Array{Float64}
    s_vec = zeros(OrderXY+1)::Array{Float64}
    sdot_vec = zeros(OrderXY+1)::Array{Float64}

    for i = 1:OrderXY+1
            s_vec[i] =s^(OrderXY-i+1)
            s0_vec[i] =s0^(OrderXY-i+1)
    end
    for i = 1:OrderXY
            sdot_vec[i] = (OrderXY+1-i)*s^(OrderXY-i)
    end

    XCurve0 = dot(coeffX,s0_vec)
    YCurve0 = dot(coeffY,s0_vec)

    XCurve  = dot(coeffX,s_vec)
    YCurve  = dot(coeffY,s_vec)

    dX = dot(coeffX,sdot_vec)
    dY = dot(coeffY,sdot_vec)      

    xyVectorAngle = atan2(y-YCurve0,x-XCurve0)
    xyPathAngle = atan2(dY,dX)

    

    ey = eyabs*sign(sin(xyVectorAngle-xyPathAngle))

    epsi = mod((psi+pi),(2*pi))-pi-xyPathAngle
    epsi = mod((epsi+pi),(2*pi))-pi

    #T Calcuate the error due to the conversion in the curvilinear abscissa
    yBack = YCurve + ey*cos(xyPathAngle)
    xBack = XCurve - ey*sin(xyPathAngle)
    Error = sqrt((y-yBack)^2 + (x-xBack)^2)
    if Error[1] >= 0.001
        warn("problem with approximation of x and y pos. Error: $Error, i: $itercount")
        # println("Distance New= ",DistanceNew)
        # println("S_Value= ",S_Value)
        # println("x= ",x)
        # println("y= ",y)
        # println("coeffX= ",coeffX)
        # println("coeffY= ",coeffY)
        # println("j_vec= ",j_vec)
        # println("s_nearest= ",s_nearest)
        # println("j= ",(s_nearest-ds):discretization:(s_nearest+ds))
    end
    # #endT

    zCurr_s = zeros(4)
    v_abs = sqrt(states_x[3].^2 + states_x[4].^2)
    zCurr_s = [s ey epsi v_abs]
    return zCurr_s, coeffCurv
end