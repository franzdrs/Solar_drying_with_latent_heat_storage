   # -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 11:35:23 2022

@author: FranzR
"""

from numpy import log, zeros, linspace, ones, where, amax, sum

#########################
# SUBROUTINE TO READ THE INPUT DATA FILE AND ASSIGN EACH DATA VARIABLE TO A LIST WHOSE LENGTH IS THE NUMBER OF DATA ROWS.
# AS IT IS, THE SUBROUTINE IS MADE TO READ THE OUPUT FILE FROM THE SAH PROGRAM WHICH HAS 6 ROWS WITH SETUP INFO, WHICH ARE SKIPPED
# THE USER CAN MODIFY THE SUBROUTINE OR THE INPUT FILE TO MATCH IF THE DATA FILE FORMAT IS DIFFERENT
def extract_data(filename):
    infile = open(filename, 'r')
    infile.readline()   #skip the first line
    infile.readline()   #skip the second line
    infile.readline()   #skip the third line
    infile.readline()   #skip the fourth line
    infile.readline()   #skip the fifth line
    infile.readline()   #skip the sixth line
    date = []
    hour_min = []
    hour_dec = []
    time_min = []
    I = []
    Tamb = []
    wind = []
    Patm = []
    X_in = []
    Ta_in = []
    mdot_a = []
        
    for line in infile:
        words = line.split()    
        
        # words[0] is date, words[1] is hour_min, words[2] is time, etc
        date.append(str(words[0]))
        hour_min.append(str(words[1]))
        hour_dec.append(float(words[2]))
        time_min.append(float(words[3]))
        I.append(float(words[4]))
        Tamb.append(float(words[5]))
        wind.append(float(words[6]))
        Patm.append(float(words[7]))
        X_in.append(float(words[8]))
        Ta_in.append(float(words[9]))
        mdot_a.append(float(words[10]))
        
    infile.close()
        
    return date, hour_min, hour_dec, time_min, I, Tamb, wind, Patm, X_in, Ta_in, mdot_a
#################################



# MAIN PROGRAM

Input_file = "SAH_Output_20m2_Stutt_7_8Aug2022 as Input.txt"
Output_file = "LHS_Output_SAH_20m2_Slab3m_Stutt_7_8Aug2022.txt"
Output_file_DETAILS_T = "LHS_Output_SAH_20m2_Slab3m_Stutt_7_8Aug2022_DETAILS_T.txt"
Output_file_DETAILS_H = "LHS_Output_SAH_20m2_Slab3m_Stutt_7_8Aug2022_DETAILS_H.txt"

# READING DATA FROM FILE
date, hour_min, hour_dec, t_min, Ins, Tamb, wind, P_atm, X_in, Ta_in, mdot_a = extract_data(Input_file)

### USER INPUTS ###

Data_time_step = 60.    # This time resolution of the data in the input file (s)
relax = 1.               # Underrelaxation factor, can be lower than 1 for stability if required, but more than 0.

# PCM PROPERTIES
rho_s = 1517.               # solid PCM density (kg/m3)
rho_l = 1442.               # liquid PCM density (kg/m3) 

ks = 1.09                  # solid PCM conductivity (W/m-K)
kl = 0.54                  # liquid PCM conductivity (W/m-K)
                 
cs = 2000.              # solid PCM specific heat (J/kg-K)
cl = 2000.              # liquid PCM specific heat (J/kg-K)

Tm = 38.0               # Midpoint of PCM melting temperature range (°C)
T_range = 4.0           # Temperature range over which the PCM melts (°C)

Hs = 0.0                # Enthalpy of solid at Ts (J/kg) (we assume that PCM enthalpy is 0 at Ts so negative below Ts) 

L = 190000.0             # PCM latent heat (J/kg)    
Hl = Hs + L	            # Enthalpy of liquid at Tl (J/kg)

# AIR PROPERTIES ASSUMED CONSTANT, THE SAME AS IN SOLAR AIR HEATER MODEL
rho_a = 1.127           # Air density (at about 40 °C) (kg/m3)
c_a = 1007.             # Air specific heat (J/kg-K)
mu_a = 1.907e-5         # Air viscosity (at 40 °C) (Pa-s)
k_a = 0.02735           # Themal conductivity of air (at 40 °C) (W/m-K)
Pr = 0.706              # Prandtl number of air

# PROPERTIES OF CONTAINER WALLS
rho_w = 2670.           # Container wall density (kg/m3)
thickn_w = 0.0005       # Container wall thickness (m)
c_w = 890.              # Container wall specific heat (J/kg-K)
k_w = 130.              # Container wall thermal conductivity (W/m-K)

# STORAGE VESSEL PROPERTIES
N_slabs = 30            # Number of slabs in storage vessel
N_channels = N_slabs    # Assuming that there is a channel between inner slabs and the outer channels adjacent to the lateral walls are half as wide. 
Slab_length = 3.0       # Length of PCM slab (m)
Slab_width_s = 0.5        # Width/height of PCM in slab when all PCM is solid (m)

Slab_thickness = 0.01   # Thickness of the slab (m)
Air_gap = 0.01          # Width of the air channels between slabs (m)

# SPACE AND TIME GRID SIZE
dt = 30.                # Timestep (s)
dx = 0.05               # Mesh size along slab length (tentative) (m)
dy = 0.001              # Mesh size along slab thickness (tentative) (m)    

# TOLERANCE FOR THE DIFFERENCE BETWEEN CONSECUTIVE ITERATIONS
TOL_delta_H = 0.1       # Tolerance of nodal enthalpy (J/kg)
TOL_delta_Ta = 0.0001   # Tolerance of nodal air temperature (°C)
TOL_delta_Tw = 0.0001   # Tolerance of nodal wall temperature (°C)


# CHOOSE HOW OFTEN TO STORE DATA (IN MINUTES)
Store_dt_min = 1.    # Data storage interval (min)
Store_dt_min_details = 5.    # Data storage interval for detailed slab data (min)

# INITIAL CONDITIONS
Tw0 = Tamb[0]           # Initial temperature of container wall (°C)
T0 = Tamb[0]            # Initial temperature of PCM (°C)
Ta0 = Tamb[0]           # Initial temperature of air inside storage (°C)

# CALCULATIONS START HERE
expansion = rho_s/rho_l
Tl = Tm + T_range/2.    # Upper bound of melting temperature range  (°C)
Ts = Tm - T_range/2.    # Lower bound of melting temperature range (°C)

Duct_width = Slab_width_s*expansion  # This would be equal to the slab width/height when the whole PCM is liquid, due to the expansion, and would
                                   # be equal to the duct width (m)
                                   
Slab_half_thickness = Slab_thickness/2.0  # The half thickness, which is the actual modeled region (m)

Mass_pcm = N_slabs*(Slab_length*Slab_width_s*Slab_thickness)*rho_s     # Total mass of PCM in LHS (kg)

Half_air_gap = Air_gap/2.       # Half thickness, which is what is taken in the model since we are considering only the heat transfer to/from one slab and using symmetry

Dh = 4*Duct_width*Air_gap/(2*(Duct_width + Air_gap))    # Hydraulic diameter of air channels, assuming the slabs go al the way to the walls (m)

# CREATING THE NECESSARY VARIABLES, VECTORS AND ARRAYS FOR THE SCHEME
M = int(round(Slab_length/dx))              # Number of divisions along slab length (x direction), -   
N = int(round(Slab_half_thickness/dy))      # Number of divisions along half the slab thickness (y direction), - 
dy = Slab_half_thickness/(N+0.5)

x_position = linspace(0, Slab_length, M+1)              # Vector of node coordinates along the slab length (m)
y_position = linspace(dy/2, Slab_half_thickness, N+1)   # Vector of node coordinates along the slab half thickness (m)

dx = x_position[1]- x_position[0]                       # Mesh size along slab lenth (actual) (m)
dy = y_position[1]- y_position[0]                       # Mesh size along slab thickness (actual) (m)

fe = (dy/2.)/(dy/2 + thickn_w)      # Ratio of distances for harmonic mean interface conductivity between container wall and PCM (Book Patankar page. 45)

Ta = zeros(M+1)                 # Vector of air temperatures along slab length (°C)
Ta_new = zeros(M+1)             # Vector of air temperatures along slab length at new time level (°C)

Ta[:] = Ta0                     # Initial air temperature in the air channel (°C)
Ta_new[:] = Ta

Tw = zeros(M+1)                 # Vector of container wall temperatures along slab length (°C)
Tw_new = zeros(M+1)             # Vector of container wall temperatures along slab length at new time level (°C)

Tw[:] = Tw0                     # Initial container wall temperature (°C)
Tw_new[:] = Tw

T = zeros((N+1,M+1))            # Matrix with temperature values in PCM grid (°C). First index is the y_direction, second index x_direction
T_new = zeros((N+1,M+1))        # Matrix with temperature values in PCM grid at the next time level (°C)

T[:] = T0                          #Initial temperature of PCM (°C)
T_new[:] = T

c = zeros((N+1,M+1))                # Matrix of specific heat values in the grid (J/kg-K)
k = zeros((N+1,M+1))                # Matrix of thermal conductivity values in the grid (W/m-K)

c = where(T >= Tm, cl,cs)           #  Initializing values of matrix c (J/kg-K)
k = where(T >= Tm, kl,ks)           #  Initializing values of matrix k (W/m-K)

H = zeros((N+1,M+1))                # Matrix with enthalpy values in PCM grid (J/kg)
H_new = zeros((N+1,M+1))            # Matrix with enthalpy values in PCM grid at the next time level (J/kg)

H = where(T <= Ts, c*(T - Ts), where(T >= Tl, L + c*(T - Tl), L*(T - Ts/(T_range))))  # Initializing values of matrix H 
H_new[:] = H

Lambda = zeros((N+1,M+1))           # Matrix to track liquid fraction, - 
Lambda = where(H <= Hs, 0.0, where(H >= Hl, 1.0, H/Hl))

grid_weights = zeros((N+1,M+1))          # Matrix with weight of each control volume to calculate true avg enthalpy and true avg lambda of slab

grid_weights[:] = 1.              # First set all elements to unit weight

j=0                               # Here we assign the corresponding weights to the control volumes.  
for i in range(0,N):
    grid_weights[i,j] = 0.5   
j=M    
for i in range(0,N):
    grid_weights[i,j] = 0.5
grid_weights[N,0] = 0.25
grid_weights[N,M] = 0.25
i=N
for j in range(1,M):
    grid_weights[i,j] = 0.5


#Lambda_avg = mean(Lambda)           # Calculates the mean liquid fraction of the entire PCM region

Lambda_avg = sum(Lambda*grid_weights)/sum(grid_weights)  # Calculates the weighted average liquid fraction of the entire PCM region
Total_Enthalpy = sum(H*grid_weights)/sum(grid_weights)*Mass_pcm   # Calculates the total enthalpy stored from the lower end of the phase change
Avg_Temperature = sum(T*grid_weights)/sum(grid_weights)  # Calculates the weighted average PCM temperature

Store_dt = Store_dt_min*60                  # Data storage interval, s
Store_dt_details = Store_dt_min_details*60  # Data storage interval for detailed slab data, s

Total_sim_t = len(t_min)*Data_time_step      # The total simulation time, s

delta_H = ones((N+1,M+1))       # Matrix of iteration changes of H in each node, inititalized to a value above tolerance desired (kJ/kg)
delta_Ta = ones(M+1)            # Matrix of iteration changes of Ta in each node, initialized to a value above tolerance desired (°C)
delta_Tw = ones(M+1)            # Matrix of iteration changes of Tw in each node, initialized to a value above tolerance desired (°C)

iterations_tot = 0
      
Curr_sim_time = 0.0                       # The current simulation time during the calculations (s)

Next_store_time = 0.                     # Initializing the variable keeping track of the next simulation time to store values
Next_store_time_details = 0.             # Initializing the variable keeping track of the next simulation time to store detail values


# HERE WE WRITE TO THE FILE THAT IS GOING TO BE USED BY A DOWNSTREAM COMPONENT AS INPUT FILE
outfileT = open(Output_file, "w")

outfileT.write( "rho_s:%-8.1f rho_l:%-8.1f ks:%-8.2f kl:%-8.2f cs:%-8.1f cl:%-8.1f L:%-8.1f Tm:%-8.1f" % (rho_s, rho_l, ks, kl, cs, cl, L, Tm))
outfileT.write("\n")
outfileT.write( "N_slabs:%-8d Slab_Length:%-8.2f Slab_width_s:%-8.2f Slab_thickness:%-8.2f PCM_mass:%-8.2f Air_gap:%-8.2f T0:%-8.1f Ta0:%-8.1f" % (N_slabs, Slab_length, Slab_width_s, Slab_thickness, Mass_pcm, Air_gap, T0, Ta0))
outfileT.write("\n")
outfileT.write( "T_range:%-8.1f Wall_thick:%-8.4f" % (T_range, thickn_w))
outfileT.write("\n")
outfileT.write( "dt:%-8.1f dx:%-8.4f dy:%-8.4f" % (dt, dx, dy))
outfileT.write("\n")

outfileT.write("%-12s" % date[0])
outfileT.write("%-9s" % hour_min[0])
outfileT.write("%-9.3f" % hour_dec[0])
outfileT.write("%-8.1f" % (Curr_sim_time/60.))
outfileT.write("%-8.1f" % Ins[0])
outfileT.write("%-8.1f" % Tamb[0])
outfileT.write("%-8.1f" % wind[0])
outfileT.write("%-8.0f" % P_atm[0])
outfileT.write("%-8.4f" % X_in[0])
outfileT.write("%-9.1f" % Ta[M])
outfileT.write("%-11.2f" % Avg_Temperature)
outfileT.write("%-13.1f" % Total_Enthalpy)
outfileT.write("%-11.3f" % Lambda_avg)
outfileT.write("%-11.8f" % mdot_a[0])
outfileT.write("\n")

outfileT.write("%-12s%-9s%-9s%-8s%-8s%-8s%-8s%-8s%-8s%-9s%-11s%-13s%-11s%-11s"  % ("Date", "Hour_min", "Hour_dec", "t_min", "Ins", "Tamb", "Wind", "Patm", "X", "Ta_out", "Avg_Tpcm", "Tot_Enthlp", "Lambda_avg", "mdot_a"))
outfileT.write("\n")
######################

# HERE WE WRITE TO THE FILES THAT ARE GOING TO CONTAIN DETAILED TEMPERATURE IN SLAB 
outfileT_details = open(Output_file_DETAILS_T, "w")

outfileT_details.write( "rho_s:%-8.1f rho_l:%-8.1f ks:%-8.2f kl:%-8.2f cs:%-8.1f cl:%-8.1f L:%-8.1f Tm:%-8.1f" % (rho_s, rho_l, ks, kl, cs, cl, L, Tm))
outfileT_details.write("\n")
outfileT_details.write( "N_slabs:%-8d Slab_Length:%-8.2f Slab_width_s:%-8.2f Slab_thickness:%-8.2f PCM_mass:%-8.2f Air_gap:%-8.2f T0:%-8.1f Ta0:%-8.1f" % (N_slabs, Slab_length, Slab_width_s, Slab_thickness, Mass_pcm, Air_gap, T0, Ta0))
outfileT_details.write("\n")
outfileT_details.write( "T_range:%-8.1f Wall_thick:%-8.4f dt:%-8.1f dx:%-8.4f dy:%-8.4f" % (T_range, thickn_w, dt, dx, dy))
outfileT_details.write("\n")

outfileT_details.write("%-10s" % "time,s:")
outfileT_details.write("%-10s" % Curr_sim_time)
outfileT_details.write("\n")

outfileT_details.write("%-10s" % "y,m/x,m")
for x in x_position:
    outfileT_details.write("%-8.3f" % x)
outfileT_details.write("\n")

outfileT_details.write("%-10s" % "Ta")
for value in Ta:
    outfileT_details.write("%-8.3f" % value)
outfileT_details.write("\n")

outfileT_details.write("%-10s" % "Tw")
for value in Tw:
    outfileT_details.write("%-8.3f" % value)
outfileT_details.write("\n")

for i in range(0,N+1):
    outfileT_details.write("%-10.5f" % y_position[i])
    for value in T[i,:]:
        outfileT_details.write("%-8.3f" % value)
    outfileT_details.write("\n")
outfileT_details.write("\n")
################################


# HERE WE WRITE TO THE FILES THAT ARE GOING TO CONTAIN DETAILED ENTHALPY IN SLAB 
outfileH_details = open(Output_file_DETAILS_H, "w")

outfileH_details.write( "rho_s:%-8.1f rho_l:%-8.1f ks:%-8.2f kl:%-8.2f cs:%-8.1f cl:%-8.1f L:%-8.1f Tm:%-8.1f" % (rho_s, rho_l, ks, kl, cs, cl, L, Tm))
outfileH_details.write("\n")
outfileH_details.write( "N_slabs:%-8d Slab_Length:%-8.2f Slab_width_s:%-8.2f Slab_thickness:%-8.2f PCM_mass:%-8.2f Air_gap:%-8.2f T0:%-8.1f Ta0:%-8.1f" % (N_slabs, Slab_length, Slab_width_s, Slab_thickness, Mass_pcm, Air_gap, T0, Ta0))
outfileH_details.write("\n")
outfileH_details.write( "T_range:%-8.1f Wall_thick:%-8.4f dt:%-8.1f dx:%-8.4f dy:%-8.4f" % (T_range, thickn_w, dt, dx, dy))
outfileH_details.write("\n")

outfileH_details.write("%-10s" % "time,s:")
outfileH_details.write("%-10s" % Curr_sim_time)
outfileH_details.write("\n")

outfileH_details.write("%-10s" % "y,m/x,m")
for x in x_position:
    outfileH_details.write("%-10.3f" % x)
outfileH_details.write("\n")

outfileH_details.write("%-10s" % "Ta")
for value in Ta:
    outfileH_details.write("%-10.3f" % value)
outfileH_details.write("\n")

outfileH_details.write("%-10s" % "Tw")
for value in Tw:
    outfileH_details.write("%-10.3f" % value)
outfileH_details.write("\n")

for i in range(0,N+1):
    outfileH_details.write("%-10.5f" % y_position[i])
    for value in H[i,:]:
        outfileH_details.write("%-10.1f" % value)
    outfileH_details.write("\n")
outfileH_details.write("\n")
#################################


A = 1/(rho_a*c_a*Half_air_gap)        # Factor used in equation for the air to simplify writing 
      
Next_store_time += Store_dt             # set next storage time
Next_store_time_details += Store_dt_details     # set next storage time in detailed output file

Count_stor_dt = 1        # This counts how many timesteps have passed before storing results in the first output file, so that 
                         # the avg of outlet Ta during the store_dt can be calculated.
                         # Starts at 1 because the least is to write data every time step.


# TIME LOOPS STARTS HERE
while Curr_sim_time < Total_sim_t:
    iterations = 0
          
    Tair_in = Ta_in[int(Curr_sim_time/Data_time_step)]      # Get current inlet air temperature
    mdot_air = mdot_a[int(Curr_sim_time/Data_time_step)]    # Get current air mass flow rate, kg/s
    V_a = mdot_air/rho_a                                    # Airflow rate (m3/s)
    va =  V_a/(N_channels*Air_gap*Duct_width)               # Mean air velocity through gaps, m/s

    Re = rho_a*va*Dh/mu_a                                   # Reynolds number in air channels 
    ff = (0.79*log(Re) - 1.64)**-2                          # Friction factor to be used in Nu calculation with Gnielinski eq.

    # CALCULATE THE HEAT TRANSFER COEFFICIENT BETWEEN AIR AND SURFACE;
    if Re <= 2300.0:
        Xaster = (Slab_length/Dh)/(Re*Pr)                      # Factor for Nu in laminar flow 
        Nu = 7.55 + (0.024*Xaster**(-1.14))/(1 + 0.0358*(Pr**0.17)*Xaster**(-0.64)) # Nu in parallel plates simultaneously developing flow ROHSENOW PAG 5.63
    else:
        Nu = (ff/8)*(Re-1000)*Pr/(1 + 12.7*(ff/8)**0.5*(Pr**(2/3)-1))    # Nu by Gnielinski (appears in practically all sources as a good eq., used by Dolado and others)

    h_aw = Nu*k_a/Dh
    
    while (amax(delta_H) > TOL_delta_H) or (amax(delta_Ta) > TOL_delta_Ta) or (amax(delta_Tw) > TOL_delta_Tw):  # This is the iteration control, the tolerance value can be varied
        
        # LOOP OVER AIR NODES FROM i=0 TO M                                
        Ta_iter = Tair_in
        Ta_new[0] = Ta_iter
        delta_Ta[0] = Ta_iter - Ta_new[0]
        for j in range (1,M+1):
            Ta_iter = (Ta[j] + (va*dt/dx)*Ta_new[j-1] + dt*A*h_aw*Tw_new[j])/(1 + dt*(va/dx + A*h_aw))
            
            delta_Ta[j] = abs(Ta_iter - Ta_new[j])
            
            Ta_new[j] = Ta_iter
                     
        # LOOP OVER CONTAINER WALL NODES
        for j in range(0,M+1):
            
            jm1 = j-1 if j > 0 else j+1         # This is to be able to use the same equations for internal and boundary nodes
            jp1 = j+1 if j < M else j-1    
            
            k_int = 1/((1-fe)/k_w + fe/(k[0,j]))  # Harmonic mean interface conducticvity for interface between container wall and PCM (Patankar). k[0,j] is the conduct of the PCM adjacent to wall 
            
            W = rho_w*thickn_w*c_w + (2*dt*k_w*thickn_w)/(dx**2) + dt*h_aw + dt*k_int/(thickn_w + dy/2) # Denominator in equation for container wall
            
            Tw_iter = (rho_w*thickn_w*c_w*Tw[j] + (dt*k_w*thickn_w/dx**2)*(Tw_new[jp1] + Tw_new[jm1]) + dt*h_aw*Ta_new[j] + (dt*k_int/(thickn_w + dy/2))*T_new[0,j])/W
            
            Tw_iter = relax*Tw_iter + (1-relax)*Tw_new[j]             # Applying underrelaxation 
            
            delta_Tw[j] = abs(Tw_iter - Tw_new[j])          # Recording the change in temperature from past iteration
            
            Tw_new[j] = Tw_iter
            
        # LOOP OVER THE NODES OF THE SURFACE ROW (i=0) FROM [0,0] UNTIL [0,M]
        i = 0
        for j in range(0,M+1):
            
            jm1 = j-1 if j > 0 else j+1         # This is to be able to use the same equations for internal and boundary nodes
            jp1 = j+1 if j < M else j-1
            
            k_int = 1/((1-fe)/k_w + fe/(k[0,j]))  # Harmonic mean interface conducticvity for interface between container wall and PCM (Patankar). k[0,j] is the conduct of the PCM adjacent to wall
           
            kx_ant = 2*k[i,jm1]*k[i,j]/(k[i,jm1]+k[i,j])      # Patankar formulas for interface conductivity using harmonic mean
            kx_post = 2*k[i,jp1]*k[i,j]/(k[i,jp1]+k[i,j])            
            ky_post = 2*k[i+1,j]*k[i,j]/(k[i+1,j]+k[i,j])                
            
            if H_new[i,j] <= Hs:
                S = rho_s/dt
            elif Hs < H_new[i,j] < Hl:
                S = ((rho_l - rho_s)/(Hl - Hs)*H_new[i,j] + rho_s)/dt
            else:
                S = rho_l/dt
            
            R = S*H[i,j] + k_int/((thickn_w + dy/2)*dy)*Tw_new[j] + (kx_post/dx**2)*T_new[i,jp1] + (kx_ant/dx**2)*T_new[i,jm1] + (ky_post/dy**2)*T_new[i+1,j]
            
            F = (k_int/((thickn_w + dy/2)*dy) + kx_post/dx**2 + kx_ant/dx**2 + ky_post/dy**2)
                       
            H_iter = (R + F*(Hs/cs - Ts))/(S + F/cs)
            H_iter = relax*H_iter + (1-relax)*H_new[i,j]             # Applying underrelaxation
            T_iter = (H_iter - Hs)/cs + Ts

            if Hs <= H_iter <= Hl:
                H_iter = (R + F*Hs*T_range/L - F*Ts)/(S + F*T_range/L)
                H_iter = relax*H_iter + (1-relax)*H_new[i,j]             # Applying underrelaxation
                T_iter = Ts + (H_iter - Hs)*T_range/L
                
            if H_iter > Hl:
                H_iter = (R + F*(Hl/cl - Tl))/(S + F/cl)
                H_iter = relax*H_iter + (1-relax)*H_new[i,j]             # Applying underrelaxation
                T_iter = (H_iter - Hl)/cl + Tl
            
            delta_H[i,j] = abs(H_iter - H_new[i,j])    # Recording the change in enthalpy from past iteration
            
            H_new[i,j] = H_iter                 # Updating value of H_new[i,j] and T_new[i,j]
            T_new[i,j] = T_iter            
                
            if H_new[i,j] <= Hs:
                k[i,j] = ks
                c[i,j] = cs
            elif Hs < H_new[i,j] < Hl:
                k[i,j] = (kl - ks)/(Hl - Hs)*H_new[i,j] + ks
                c[i,j] = (cl - cs)/(Hl - Hs)*H_new[i,j] + cs
            else:
                k[i,j] = kl
                c[i,j] = cl
               
        #NOW WE MOVE TO THE INNER ROWS, FROM ROW i=1 TO i=N-1, AND WITHIN EACH ROW MOVING FROM COLUMN j=0 TO j=M 
        for i in range(1,N):
            for j in range(0,M+1):            
                
                jm1 = j-1 if j > 0 else j+1         # This is to be able to use the same equations for internal and boundary nodes
                jp1 = j+1 if j < M else j-1
                
                kx_ant = 2*k[i,jm1]*k[i,j]/(k[i,jm1]+k[i,j])      # Patankar formulas for interface conductivity using harmonic mean
                kx_post = 2*k[i,jp1]*k[i,j]/(k[i,jp1]+k[i,j])
                ky_ant = 2*k[i-1,j]*k[i,j]/(k[i-1,j]+k[i,j])
                ky_post = 2*k[i+1,j]*k[i,j]/(k[i+1,j]+k[i,j])
                
                if H_new[i,j] <= Hs:
                    S = rho_s/dt
                elif Hs < H_new[i,j] < Hl:
                    S = ((rho_l - rho_s)/(Hl - Hs)*H_new[i,j] + rho_s)/dt
                else:
                    S = rho_l/dt

                R = S*H[i,j] + (kx_post/dx**2)*T_new[i,jp1] + (kx_ant/dx**2)*T_new[i,jm1] + (ky_post/dy**2)*T_new[i+1,j] + (ky_ant/dy**2)*T_new[i-1,j]
                
                F = (kx_post/dx**2 + kx_ant/dx**2 + ky_post/dy**2 + ky_ant/dy**2)
                
                H_iter = (R + F*(Hs/cs - Ts))/(S + F/cs)
                H_iter = relax*H_iter + (1-relax)*H_new[i,j]             # Applying underrelaxation
                T_iter = (H_iter - Hs)/cs + Ts
                
                if Hs <= H_iter <= Hl:
                    H_iter = (R + F*Hs*T_range/L - F*Ts)/(S + F*T_range/L) 
                    H_iter = relax*H_iter + (1-relax)*H_new[i,j]             # Applying underrelaxation
                    T_iter = Ts + (H_iter - Hs)*T_range/L
                                   
                if H_iter > Hl:
                    H_iter = (R + F*(Hl/cl - Tl))/(S + F/cl)
                    H_iter = relax*H_iter + (1-relax)*H_new[i,j]             # Applying underrelaxation
                    T_iter = (H_iter - Hl)/cl + Tl
                    
                delta_H[i,j] = abs(H_iter - H_new[i,j])         # Recording the change in enthalpy from past iteration
        
                H_new[i,j] = H_iter     # Updating value of H_new[i,j] and T_new[i,j]
                T_new[i,j] = T_iter
                                
                if H_new[i,j] <= Hs:
                    k[i,j] = ks
                    c[i,j] = cs
                elif Hs < H_new[i,j] < Hl:
                    k[i,j] = (kl - ks)/(Hl - Hs)*H_new[i,j] + ks
                    c[i,j] = (cl - cs)/(Hl - Hs)*H_new[i,j] + cs
                else:
                    k[i,j] = kl
                    c[i,j] = cl
        
        #NOW THE SYMMETRY ROW i=N, GOING FROM COLUMN J=0 TO J
        i = N
        for j in range(0,M+1):
           
            jm1 = j-1 if j > 0 else j+1         # This is to be able to use the same equations for internal and boundary nodes
            jp1 = j+1 if j < M else j-1                              
            
            kx_ant = 2*k[i,jm1]*k[i,j]/(k[i,jm1]+k[i,j])      # Patankar formulas for interface conductivity using harmonic mean
            kx_post = 2*k[i,jp1]*k[i,j]/(k[i,jp1]+k[i,j])
            ky_ant = 2*k[i-1,j]*k[i,j]/(k[i-1,j]+k[i,j])

            if H_new[i,j] <= Hs:
                S = rho_s/dt
            elif Hs < H_new[i,j] < Hl:
                S = ((rho_l - rho_s)/(Hl - Hs)*H_new[i,j] + rho_s)/dt
            else:
                S = rho_l/dt     
    
            R = S*H[i,j] + (kx_post/dx**2)*T_new[i,jp1] + (kx_ant/dx**2)*T_new[i,jm1] + 2*(ky_ant/dy**2)*T_new[i-1,j]    
    
            F = (kx_post/dx**2 + kx_ant/dx**2 + 2*ky_ant/dy**2)
        
            H_iter = (R + F*(Hs/cs - Ts))/(S + F/cs)
            H_iter = relax*H_iter + (1-relax)*H_new[i,j]             # Applying underrelaxation
            T_iter = (H_iter - Hs)/cs + Ts
            
            if Hs <= H_iter <= Hl:
                H_iter = (R + F*Hs*T_range/L - F*Ts)/(S + F*T_range/L) 
                H_iter = relax*H_iter + (1-relax)*H_new[i,j]             # Applying underrelaxation
                T_iter = Ts + (H_iter - Hs)*T_range/L
            
            if H_iter > Hl:
                H_iter= (R + F*(Hl/cl - Tl))/(S + F/cl)
                H_iter = relax*H_iter + (1-relax)*H_new[i,j]             # Applying underrelaxation
                T_iter= (H_iter - Hl)/cl + Tl
                
            delta_H[i,j] = abs(H_iter - H_new[i,j])                 # Recording the change in enthalpy from past iteration
            
            H_new[i,j] = H_iter         # Updating value of H_new[i,j] and T_new[i,j]
            T_new[i,j] = T_iter
                
            if H_new[i,j] <= Hs:
                k[i,j] = ks
                c[i,j] = cs
            elif Hs < H_new[i,j] < Hl:
                k[i,j] = (kl - ks)/(Hl - Hs)*H_new[i,j] + ks
                c[i,j] = (cl - cs)/(Hl - Hs)*H_new[i,j] + cs
            else:
                k[i,j] = kl
                c[i,j] = cl
                        
        iterations += 1
        iterations_tot += 1
  
    
    print ("max_delta_Ta", amax(delta_Ta))
    print ("max_delta_Tw", amax(delta_Tw))
    print ("max_delta_H", amax(delta_H))
    
    print ("t", Curr_sim_time)
    print ("iterations", iterations)       
    
    delta_H[:] = 1.0              # Reset the values of delta_H and delta_Ta in each node to a value above the tolerance desired
    delta_Ta[:] = 1.0
    delta_Tw[:] = 1.0
    
     
    Lambda = where(H_new <= Hs, 0.0, where(H_new >= Hl, 1.0, H_new/Hl))
    #Lambda_avg = mean(Lambda)
    pcm_width = Slab_width_s + Lambda_avg*(Duct_width - Slab_width_s)                   # This PCM width is the actual width/height of the pcm inside the slabs as it melts
                                                                                    # according to how much pcm has melted (liquid fraction). it can go from the level when
                                                                                    # all pcm is solid, to the level when all pcm is liquid, according to the calculated 
                                                                                    # average liquid fraction
    
    Ta_out_mixed = (Ta_new[M]*pcm_width + Tair_in*(Duct_width - pcm_width))/Duct_width      # This calculates the mixed temperature of the air at the outlet of the slab and the 
                                                                                          # air which was above the pcm level, which is leaves de LHS as it arrived, with the inlet T
                                                                                          # which is the collector outlet temperature.
      
    Curr_sim_time += dt             # Advance time step
        
    H[:] = H_new                    # Replacing the new calculated enthalpy and temperatures into the arrays of current values to move into next timestep
    T[:] = T_new
    Ta[:] = Ta_new
    Tw[:] = Tw_new
    
    Lambda_avg = sum(Lambda*grid_weights)/sum(grid_weights)  # Calculates the weighted average liquid fraction of the entire PCM region
    Total_Enthalpy = sum(H*grid_weights)/sum(grid_weights)*Mass_pcm   # Calculates the total enthalpy stored from the lower end of the phase change
    Avg_Temperature = sum(T*grid_weights)/sum(grid_weights)  # Calculates the weighted average PCM temperature
    
    # STORING DATA IN OUTPUT FILES
 
    index = int(Curr_sim_time/Data_time_step)
    if Curr_sim_time % Data_time_step == 0:
       index -= 1
       
    if (Next_store_time - Curr_sim_time) < dt:
        outfileT.write("%-12s" % date[index])
        outfileT.write("%-9s" % hour_min[index])
        outfileT.write("%-9.3f" % hour_dec[index])
        outfileT.write("%-8.1f" % (Curr_sim_time/60.))
        outfileT.write("%-8.1f" % Ins[index])
        outfileT.write("%-8.1f" % Tamb[index])
        outfileT.write("%-8.1f" % wind[index])
        outfileT.write("%-8.0f" % P_atm[index])
        outfileT.write("%-8.4f" % X_in[index])
        outfileT.write("%-9.3f" % Ta_out_mixed)
        outfileT.write("%-11.2f" % Avg_Temperature)
        outfileT.write("%-13.1f" % Total_Enthalpy)
        outfileT.write("%-11.3f" % Lambda_avg)
        outfileT.write("%-11.8f" % mdot_a[index])
        outfileT.write("\n")
        
        print ("saving....")
    
        Count_stor_dt = 0 
        Next_store_time += Store_dt
    
    if (Next_store_time_details - Curr_sim_time) < dt:

        outfileT_details.write("%-10s" % "time,s:")
        outfileT_details.write("%-10s" % Curr_sim_time)
        outfileT_details.write("\n")
        
        outfileT_details.write("%-10s" % "y,m/x,m")
        for x in x_position:
            outfileT_details.write("%-8.3f" % x)
        outfileT_details.write("\n")

        outfileT_details.write("%-10s" % "Ta")
        for value in Ta_new:
            outfileT_details.write("%-8.3f" % value)
        outfileT_details.write("\n")
        
        outfileT_details.write("%-10s" % "Tw")
        for value in Tw_new:
            outfileT_details.write("%-8.3f" % value)
        outfileT_details.write("\n")

        for i in range(0,N+1):
            outfileT_details.write("%-10.5f" % y_position[i])
            for value in T[i,:]:
                outfileT_details.write("%-8.3f" % value)
            outfileT_details.write("\n")
        outfileT_details.write("\n")

        outfileH_details.write("%-10s" % "time,s:")
        outfileH_details.write("%-10s" % Curr_sim_time)
        outfileH_details.write("\n")
        
        outfileH_details.write("%-10s" % "y,m/x,m")
        for x in x_position:
            outfileH_details.write("%-10.3f" % x)
        outfileH_details.write("\n")

        outfileH_details.write("%-10s" % "Ta")
        for value in Ta_new:
            outfileH_details.write("%-10.3f" % value)
        outfileH_details.write("\n")
        
        outfileH_details.write("%-10s" % "Tw")
        for value in Tw_new:
            outfileH_details.write("%-11.3f" % value)
        outfileH_details.write("\n")

        for i in range(0,N+1):
            outfileH_details.write("%-10.5f" % y_position[i])
            for value in H[i,:]:
                outfileH_details.write("%-10.1f" % value)
            outfileH_details.write("\n")
        outfileH_details.write("\n")
             
        Next_store_time_details += Store_dt_details
        
    Count_stor_dt += 1

outfileT_details.close()
outfileH_details.close()
outfileT.close()


