function trackFrameConversion(states_x::Array{Float64},x_track::Array{Float64},y_track::Array{Float64},trackCoeff::classes.TrackCoeff, itercount::Int64, N::Int64, dt::Float64 )

    OrderXY        = trackCoeff.nPolyXY
    OrderThetaCurv = trackCoeff.nPolyCurvature
    ds             = trackCoeff.ds

    # grab current states of the vehicle
    x       = states_x[1]
    y       = states_x[2]
    psi     = states_x[5]
    v_abs   = sqrt(states_x[3].^2 + states_x[4].^2)

    nodes_center = [x_track; y_track]

    N_nodes_poly_back  = convert(Int64,ceil(0.5*v_abs*N*dt/trackCoeff.ds))
    N_nodes_poly_front = convert(Int64,ceil(1.5*v_abs*N*dt/trackCoeff.ds))
    nPoints            = N_nodes_poly_back + N_nodes_poly_front

    dist_vec = sqrt((x-x_track).^2+(y-y_track).^2)

    dist_lane, idx_min = findmin(dist_vec)

    N_nodes_center  = size(nodes_center,2) # how many nodes define the track 
    ind_start       = idx_min-N_nodes_poly_back 
    ind_end         = idx_min+N_nodes_poly_front 

    if idx_min <= N_nodes_poly_back 
        nodes_near = hcat(nodes_center[:,N_nodes_center+ind_start:N_nodes_center],nodes_center[:,1:ind_end]) 
    elseif idx_min+N_nodes_poly_front >= N_nodes_center
        nodes_near = hcat(nodes_center[:,ind_start:N_nodes_center],nodes_center[:,1:ind_end-N_nodes_center]) 
    else 
        nodes_near = nodes_center[:,ind_start:ind_end]
    end

    s_interp_start = (ind_start-1)*ds 
    s_nearest      = (idx_min-1)*ds
    
    nodes_near_X = vec(nodes_near[1,:])
    nodes_near_Y = vec(nodes_near[2,:])

    index = 0 

    itpMatrix = zeros(nPoints+1,OrderXY+1)
   

    for i = s_interp_start:ds:s_interp_start+(nPoints*ds)
        for k = 0:OrderXY
            index = convert(Int,(i-s_interp_start)/ds+1)
            itpMatrix[index,OrderXY+1-k] = i^k
        end
    end

    Matrix4th = zeros(nPoints+1,OrderThetaCurv+1) #generate a matrix of 4th order to approximate Theta and the curvature
    for i = s_interp_start:ds:s_interp_start+nPoints*ds
        for k = 0:OrderThetaCurv
            index = convert(Int,(i-s_interp_start)/ds+1)
            Matrix4th[index,OrderThetaCurv+1-k] = i^k
        end
    end

    coeffY = itpMatrix\nodes_near_Y
    coeffX = itpMatrix\nodes_near_X

    Counter        = 1
    Discretization = 0.01*ds
    j_vec          = zeros(OrderXY+1) 


    if s_nearest >= ds
        S_Value = zeros(0:Discretization:2*ds) 
        DistanceNew = zeros(0:Discretization:2*ds)
        for j=(s_nearest-ds):Discretization:(s_nearest+ds) 
            
            for i = 1:OrderXY+1
                j_vec[i] =j^(OrderXY+1-i)
            end
            XCurve = dot(coeffX, j_vec)
            YCurve = dot(coeffY, j_vec)

            S_Value[Counter] = j
            DistanceNew[Counter] = sqrt((x-XCurve).^2+(y-YCurve).^2) 
            Counter = Counter + 1
        end
    else 
        S_Value = zeros(0:Discretization:1*ds) 
        DistanceNew = zeros(0:Discretization:1*ds)
        for j=s_nearest:Discretization:(s_nearest+ds) 
    
            for i = 1:OrderXY+1
                j_vec[i] =j^(OrderXY+1-i)
            end
            XCurve = dot(coeffX, j_vec)
            YCurve = dot(coeffY, j_vec)

            S_Value[Counter] = j
            DistanceNew[Counter] = sqrt((x-XCurve).^2+(y-YCurve).^2) 
            Counter = Counter + 1
        end
    end

    eyabs, idx_min_Dist = findmin(DistanceNew)

    s = S_Value[idx_min_Dist]

    if s >= 0.01*ds
        s0 = s-0.01*ds
        s0_vec = zeros(OrderXY+1,1)
        s_vec = zeros(OrderXY+1,1)
        sdot_vec = zeros(OrderXY+1,1)::Array{Float64,2}

        for i = 1:OrderXY+1
                s_vec[i] =s^(OrderXY+1-i)
                s0_vec[i] =s0^(OrderXY+1-i)
                sdot_vec[i] = (OrderXY+1-i)*s^(OrderXY-i)
        end

        XCurve0 = coeffX'*s0_vec
        YCurve0 = coeffY'*s0_vec

        XCurve = coeffX'*s_vec
        YCurve = coeffY'*s_vec
        dX = coeffX'*sdot_vec 
        dY = coeffY'*sdot_vec      

        xyVectorAngle = atan2(y-YCurve0,x-XCurve0)
        xyPathAngle = atan2(dY,dX)
    else
        
        s0 = s+0.01*ds
        s0_vec = zeros(OrderXY+1,1)
        s_vec = zeros(OrderXY+1,1)
        sdot_vec = zeros(OrderXY+1,1)::Array{Float64,2}

        for i = 1:OrderXY+1
                s_vec[i] =s^(OrderXY+1-i)
                s0_vec[i] =s0^(OrderXY+1-i)
               
        end
        sdot_vec = [6*s^5 ;5*s^4; 4*s^3; 3*s^2; 2*s ;1; 0]

        XCurve0 = coeffX'*s0_vec
        YCurve0 = coeffY'*s0_vec

        XCurve = coeffX'*s_vec
        YCurve = coeffY'*s_vec
        
        dX = dot(coeffX,sdot_vec) 
        dY = dot(coeffY,sdot_vec)

        xyVectorAngle = atan2(y-YCurve0,x-XCurve0)
        xyPathAngle = atan2(dY,dX)

    end

    ey = eyabs*sign(sin(xyVectorAngle-xyPathAngle))


    epsi = mod((psi+pi),(2*pi))-pi-xyPathAngle
    epsi = mod((epsi+pi),(2*pi))-pi



   
    yBack = YCurve + ey*cos(xyPathAngle)
    xBack = XCurve - ey*sin(xyPathAngle)
    Error = sqrt((y-yBack)^2 + (x-xBack)^2)
    if Error[1] >= 0.001
        warn("problem with approximation of x and y pos. Error: $Error, i: $itercount")
    end


    b_curvature_vector = zeros(nPoints+1)

    Counter = 1


    for j = s_interp_start:ds:s_interp_start+nPoints*ds
        dX = dot(coeffX,[6*j^5, 5*j^4, 4*j^3, 3*j^2, 2*j, 1, 0]) 
        dY = dot(coeffY,[6*j^5, 5*j^4, 4*j^3, 3*j^2, 2*j, 1, 0])
        ddX = dot(coeffX,[30*j^4, 20*j^3, 12*j^2, 6*j, 2, 0, 0])
        ddY = dot(coeffY,[30*j^4, 20*j^3, 12*j^2, 6*j, 2, 0, 0])
                       
        curvature = (dX*ddY-dY*ddX)/(dX^2+dY^2)^(3/2) 
        b_curvature_vector[Counter] = curvature 
           
        Counter = Counter + 1
    end
  

    coeffCurv  = Matrix4th\b_curvature_vector

        j = s
        dX = dot(coeffX,[6*j^5, 5*j^4, 4*j^3, 3*j^2, 2*j, 1, 0]) 
        dY = dot(coeffY,[6*j^5, 5*j^4, 4*j^3, 3*j^2, 2*j, 1, 0])
        ddX = dot(coeffX,[30*j^4, 20*j^3, 12*j^2, 6*j, 2, 0, 0])
        ddY = dot(coeffY,[30*j^4, 20*j^3, 12*j^2, 6*j, 2, 0, 0])

  
    zCurr_s = zeros(4)
    zCurr_s = [s ey epsi v_abs]
    return zCurr_s, coeffCurv
end

