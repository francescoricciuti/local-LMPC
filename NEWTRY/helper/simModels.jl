# This file contains the models used to simulate the vehicle



function simModel_exact_dyn_x(z::Array{Float64},u::Array{Float64},modelParams::classes.ModelParams) # this function simulates an exact model (no model mismatch)
    # define usefull parameters
    dt = modelParams.dt    # time step
    dtn = dt/10            # simulation time step
    t = 0:dtn:dt           # vector of times for the simulation
    z_final = copy(z)      # initialize z_final

    for i=1:length(t)-1
        z_final = simModel_dyn_x(z_final,u,dtn,modelParams,i) 
    end
    return z_final         # this will be the output of the function (the one that in "main" is called "zCurr_x[i+1,:]")
end


function simModel_dyn_x(z::Array{Float64},u::Array{Float64},dt::Float64,modelParams::classes.ModelParams,exact_sim_i::Int64)


    zNext::Array{Float64}     # initialize the z_next ...
    zNext   = z               # ... and set it equal to 'z', which is the state in which we are when this simulator is called ("zCurr_x[i,:]")
    x       = z[1]            # initialize all the states with the values of 'z'
    y       = z[2]
    v_x     = z[3]
    v_y     = z[4]
    psi     = z[5]
    psi_dot = z[6]

    # Initialize model parameters
  
    l_A     = modelParams.l_A
    l_B     = modelParams.l_B
    max_alpha = modelParams.max_alpha
    m         = modelParams.mass
    mu        = modelParams.mu
    g         = modelParams.g
    I_z       = modelParams.I_z
    B         = modelParams.B
    C         = modelParams.C

    # Determine quantities of the model

    F_xr      = m*u[1]    
    FMax      = mu*m*g / 2.0    
    if F_xr > FMax      
        F_xr = FMax    
    elseif F_xr < -FMax  
        F_xr = -FMax 
    end

    # determine slip angles      
    if v_x < 0.1        
        alpha_f = 0.0        
        alpha_r = 0.0      
    else        
        alpha_f = atan( (v_y+l_A*psi_dot) / v_x ) - u[2]        
        alpha_r = atan( (v_y-l_B*psi_dot) / v_x)       
    end
    if exact_sim_i==1 && max(abs(alpha_f),abs(alpha_r))>max_alpha/180*pi
        warn("Large slip angles: alpha_f = $(alpha_f*180/pi)°, alpha_r = $(alpha_r*180/pi)° , x =$x, y = $y")
    end
    
    F_yf = -FMax * sin(C*atan(B*alpha_f))
    F_yr = -FMax * sin(C*atan(B*alpha_r))
    
    if F_yr > sqrt(FMax^2 - F_xr^2)        
        F_yr = sqrt(FMax^2 - F_xr^2)    
    elseif  F_yr < -sqrt(FMax^2 - F_xr^2)  
        F_yr = -sqrt(FMax^2 - F_xr^2) 
    end

    # Solve one step of the model

    zNext[1] = x + dt*(v_x*cos(psi) - v_y*sin(psi))
    zNext[2] = y + dt*(v_x*sin(psi) + v_y*cos(psi))
    zNext[3] = v_x + dt*(psi_dot*v_y+1/m*(F_xr-F_yf*sin(u[2])))
    zNext[4] = v_y + dt*(-psi_dot*v_x+1/m*(F_yf*cos(u[2])+F_yr))
    zNext[5] = psi + dt*psi_dot
    zNext[6] = psi_dot + dt*(1/I_z*(l_A*F_yf*cos(u[2])-l_B*F_yr))

    # and pass the solution to "simModel_exact_dyn_x"
    return zNext
end