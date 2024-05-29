# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 11:35:23 2022

@author: FranzR
"""

from numpy import pi, exp, log, full, zeros, ones, linspace, maximum, sin, cos, amax

#########################
# SUBROUTINE TO READ THE WEATHER DATA FILE AND ASSIGN EACH DATA VARIABLE TO A LIST WHOSE LENGTH IS THE NUMBER OF DATA ROWS.
# THE USER CAN PASS THE START ROW AND END ROW TO BBE READ FROM THE DATA FILE
# THE USER CAN MODIFY THE SUBROUTINE OR THE INPUT FILE TO MATCH IF THE INPUT DATA FILE FORMAT IS DIFFERENT
def extract_data(filename, start_row, end_row):
    infile = open(filename, 'r')
    
    date = []
    hour_min = []
    time_decimal = []
    time_row = []
    I = []
    Ta = []
    rh = []
    wind = []
    Patm = []
    X = []
    Tdp = []
    mdot_a = []
        
    for row_number, line in enumerate(infile, start=1):
        if start_row <= row_number <= end_row:
        
            words = line.split()
            date.append(str(words[0]))
            hour_min.append(str(words[1]))
            time_decimal.append(float(words[2]))
            time_row.append(int(words[3]))
            I.append(float(words[4]))
            Ta.append(float(words[5]))
            rh.append(float(words[6])) 
            wind.append(float(words[7]))
            Patm.append(float(words[8]))
            mdot_a.append(float(words[9]))
            

    infile.close()
    
    for i in range(0, len(Ta)):
        Pvsat = saturation_vapor_pressure(Ta[i])
        Pv = (rh[i]/100.)*Pvsat
        X.append(0.621945*Pv/(Patm[i] - Pv))
        
        Tdp.append(6.54 + 14.526*log(Pv/1000) + 0.7389*log(Pv/1000)**2 + 0.09486*log(Pv/1000)**3 + 0.4569*(Pv/1000)**0.1984) #from ASHRAE (2017) Fundamentals pyschrometrics, eq 37
        

    return (date, hour_min, time_decimal, time_row, I, Ta, rh, wind, Patm, X, Tdp, mdot_a)
########################

# FUNCTION TO CALCULATE THE SATURATION VAPOR PRESSURE AT A GIVEN TEMPERATURE
def saturation_vapor_pressure(Tad):      # Saturation vapor pressure at some temperature, Pa (ASHRAE, 2017)
    Tabs = Tad + 273.15     
    log_Pvsat = -5.8002206e3/Tabs + 1.3914493e0 - 4.8640239e-2*Tabs + 4.1764768e-5*Tabs**2 - 1.4452093e-8*Tabs**3 + 6.5459673e0*log(Tabs)   
    Pvsat = exp(log_Pvsat)
    return Pvsat
######################


# MAIN PROGRAM
sigma = 5.6697e-8           # Stefan-Bolzmann constant (W/m2-K4)
relax = 1.0                 # Underrelaxation factor, can be lower than 1 for stability if required, but more than 0
Data_time_step = 60.        # This is the time resolution of the input data (s)



Input_file = "Data Stuttgart Aug2022 minute_LOCALTIME_mdot.txt"
Output_file = "Output_SAH_Stutt_Aug2022_27.05.2024.txt"
Output_file_DETAILS = "Output_SAH_Stutt_Aug2022_DETAILS_27.05.2024.txt"

# READING DATA FROM FILE
date, hour_min, hour_decimal, time_row, Ins, Tamb, rh_amb, wind, Patm, X, Tdewp, mdot_a = extract_data(Input_file, 12002, 14281)


### USER INPUTS ###

# COLLECTOR DIMENSIONS AND PROPERTIES
L = 6.              # Collector length  (m) 
W = 2.              # Collector width, (m)
H = 0.05            # Gap of air channel (m)
l = 0.05            # Gap between plate and cover (m) 
tilt = 25*(pi/180)   # Collector tilt (rad) (change only the first number which is the angle in degrees)

thick_g = 0.003; k_g = 1.05; c_g = 660.; rho_g = 2500.0   # Cover thickness (m), conductivity (W/m-K), specific heat (J/kg-K) and density (kg/m3)
thick_p = 0.001; k_p = 220.0; c_p = 890.; rho_p = 2700.0   # Absorber thickness (m), conductivity (W/m-K), specific heat (J/kg-K) and density (kg/m3)
thick_b = 0.001; k_b = 220.0; c_b = 890.; rho_b = 2700.0   # Back plate thickness (m), conductivity (W/m-K), specific heat (J/kg-K) and density (kg/m3)

thick_i = 0.1; k_i = 0.04                                # Back insulation thickness (m) and conductivity (J/kg-K)     
alpha_g = 0.06; tau_g = 0.84; e_g = 0.90                    # Absorptance, transmittance and emittance of cover (-)
alpha_p = 0.95; e_p_up = 0.95; e_p_down = 0.95              # Absorptance and emittance of absorber (-)
alpha_b = 0.95; e_b = 0.95                                  # Absorptance and emittance of back plate (-)                              
                                             
# AIR PROPERTIES, ASSUMED CONSTANT
rho_f = 1.127               # Air density (at about 40 °C) (kg/m3)
c_f = 1007.                 # Air specific heat (J/kg-K)
mu_f = 1.907e-5             # Air viscosity (at 40 °C) (Pa-s)
k_f = 0.02735               # Themal conductivity of air (at 40 °C) (W/m-K)
Pr = 0.706                  # Prandtl number of air

# DEFINE MESH SIZE
dx = 0.1                # Mesh size (m)
dt = 10.                # delta t desired (s)

# TOLERANCE FOR THE DIFFERENCE BETWEEN CONSECUTIVE ITERATIONS
TOL_delta_Tg = 0.00001          # Tolerance of nodal cover temperature (°C)
TOL_delta_Tp = 0.00001          # Tolerance of nodal absorber temperature (°C)
TOL_delta_Tf = 0.00001          # Tolerance of nodal fluid temperature (°C)
TOL_delta_Tb = 0.00001          # Tolerance of nodal bottom plate temperature (°C)

# CHOOSE HOW OFTEN TO STORE DATA
Store_dt_min = 1.    # Data storage interval in minutes (min)
Store_dt_min_details = 5.    # Data storage interval in minutes for detailed collector data (min)


# CALCULATIONS START HERE
     
Dh = 4*W*H/(2*(W+H))        # Hydraulic diameter of air channel (m)
Ac = L*W                    # Collector area (m2)
Ub = k_i/thick_i            # Back loss coefficient (W/m2-K)
th_diff_f = k_f/(c_f*rho_f) # Air thermal diffussivity (m2/s)

Nx = int(round(L/dx))       # Number of space intervals, -
x = linspace(0, L, Nx+1)    # Vector with the position of all nodes
dx = x[1]-x[0]              # Actual dx (m)

# Calculate factors for equations
Ag = k_g*thick_g/(dx**2); Bg = rho_g*thick_g*c_g
Af = 1/(rho_f*c_f*H)
Ap = k_p*thick_p/(dx**2); Bp = rho_p*thick_p*c_p
Ab = k_b*thick_b/(dx**2); Bb = rho_b*thick_b*c_b
 
hr_pg = zeros(Nx+1)             # Initializing vector of radiation heat transfer coefficients between absorber and cover, which are dependent on Tp and Tg
hr_gs = zeros(Nx+1)             # Initializing vector of sky radiation heat transfer coefficients which are dependent on Tg and Ts
hr_pb = zeros(Nx+1)             # Initializing vector of radiation heat transfer coefficients betweeen absorber and back plate, which are dependenet on Tp and Tb

Ra = zeros(Nx+1)
Nu_nc = zeros(Nx+1)
hnc_pg = zeros(Nx+1)

Store_dt = Store_dt_min*60  # Data storage interval in seconds (s)
Store_dt_details = Store_dt_min_details*60  # Data storage interval in seconds for detailed collector data (s)

Total_sim_time = len(time_row)*Data_time_step      # The total simulation time (s)

# Initialize main temperature vectors

Tg = full(Nx+1,Tamb[0])                 # Cover Temperature vector in current time level 
Tg_new = full(Nx+1,Tamb[0])             # Cover Temperature vector in next time level
Tp = full(Nx+1,Tamb[0] + 0.1)           # Absorber Temperature vector in current time level 
Tp_new = full(Nx+1,Tamb[0] + 0.1)       # Absorber Temperature vector in next time level
Tf = full(Nx+1,Tamb[0])                 # Fluid Temperature vector in current time level 
Tf_new = full(Nx+1,Tamb[0])             # Fluid Temperature vector in next time level
Tb = full(Nx+1,Tamb[0])                 # Back plate Temperature vector in current time level 
Tb_new = full(Nx+1,Tamb[0])             # Back plate Temperature vector in next time level

UsefulHeatTotal = 0.        # variable to sum up of useful heat over ambient T (J)
Total_I = 0.                # variable to sum up total solar radiation (J)
       
Curr_sim_time = 0.0                       # The current simulation time during the calculations, s

Next_store_time = 0.             # Initializing the variable keeping track of the next simulation time to store values
Next_store_time_details = 0.     # Initializing the variable keeping track of the next simulation time to store detail values

# OPEN FILE TO STORE OUTPUT DATA AND WRITE SETUP. THIS IS THE FILE TO BE USED BY AN DOWNSTREAM COMPONENT AS INPUT FILE
outfile1 = open(Output_file, "w")


outfile1.write( "L:%-8.3f W:%-8.3f H:%-8.3f l:%-8.3f tilt:%-8.2f thick_g:%-8.3f thick_p:%-8.3f thick_b:%-8.3f k_g:%-8.3f k_p:%-8.3f k_b:%-8.3f  " % (L, W, H, l, tilt, thick_g, thick_p, thick_b, k_g, k_p, k_b))
outfile1.write("\n")
outfile1.write( "c_g:%-8.1f c_p:%-8.1f c_b:%-8.1f rho_g:%-8.1f rho_p:%-8.1f rho_b:%-8.1f thick_i:%-8.3f k_i%-8.3f " % (c_g, c_p, c_b, rho_g, rho_p, rho_b, thick_i, k_i))
outfile1.write("\n")
outfile1.write( "alpha_g:%-8.2f alpha_p:%-8.2f alpha_b:%-8.2f tau_g:%-8.2f e_g:%-8.2f e_p_up:%-8.2f e_p_down:%-8.2f e_b:%-8.2f" % (alpha_g, alpha_p, alpha_b, tau_g, e_g, e_p_up, e_p_down, e_b))
outfile1.write("\n")
outfile1.write( "dx:%-8.2f dt:%-8.1f" % (dx, dt))
outfile1.write("\n")

outfile1.write("%-12s" % date[0])
outfile1.write("%-9s" % hour_min[0])
outfile1.write("%-9.3f" % hour_decimal[0])
outfile1.write("%-8.1f" % (Curr_sim_time/60.))
outfile1.write("%-8.1f" % Ins[0])
outfile1.write("%-8.1f" % Tamb[0])
outfile1.write("%-8.1f" % wind[0])
outfile1.write("%-8.0f" % Patm[0])
outfile1.write("%-8.4f" % X[0])
outfile1.write("%-9.1f" % Tf[Nx])
outfile1.write("%-11.8f" % mdot_a[0])
outfile1.write("\n")

outfile1.write("%-12s%-9s%-9s%-8s%-8s%-8s%-8s%-8s%-8s%-9s%-11s"  % ("Date", "Hour_min", "Hour_dec", "t_min", "Ins", "Tamb", "Wind", "Patm", "X", "Tout_avg", "mdot_a" ))
outfile1.write("\n")
##########################

# OPEN FILE TO STORE DETAILED OUTPUT DATA (COLLECTOR COMPONENTS TEMPERATURES)
outfile2 = open(Output_file_DETAILS, "w")

outfile2.write( "L:%-8.3f W:%-8.3f H:%-8.3f l:%-8.3f tilt:%-8.2f thick_g:%-8.3f thick_p:%-8.3f thick_b:%-8.3f k_g:%-8.3f k_p:%-8.3f k_b:%-8.3f  " % (L, W, H, l, tilt, thick_g, thick_p, thick_b, k_g, k_p, k_b))
outfile2.write("\n")
outfile2.write( "c_g:%-8.1f c_p:%-8.1f c_b:%-8.1f rho_g:%-8.1f rho_p:%-8.1f rho_b:%-8.1f thick_i:%-8.3f k_i%-8.3f " % (c_g, c_p, c_b, rho_g, rho_p, rho_b, thick_i, k_i))
outfile2.write("\n")
outfile2.write( "alpha_g:%-8.2f alpha_p:%-8.2f alpha_b:%-8.2f tau_g:%-8.2f e_g:%-8.2f e_p_up:%-8.2f e_p_down:%-8.2f e_b:%-8.2f" % (alpha_g, alpha_p, alpha_b, tau_g, e_g, e_p_up, e_p_down, e_b))
outfile2.write("\n")
outfile2.write( "dx:%-8.2f dt:%-8.1f" % (dx, dt))
outfile2.write("\n")

outfile2.write("%-12s%-9s%-9s%-8s%-8s%-8s%-8s"  % ("Date", "Hour_min", "Hour_dec", "t_min", "Ins", "Tamb", "Wind"))

for position in x:
    outfile2.write("G%-8.2f" % position)
for position in x:
    outfile2.write("P%-8.2f" % position)
for position in x:
    outfile2.write("F%-8.2f" % position)
for position in x:
    outfile2.write("B%-8.2f" % position)
outfile2.write("%-8s"  % "Eff")
outfile2.write("\n")    

outfile2.write("%-12s" % date[0])
outfile2.write("%-9s" % hour_min[0])
outfile2.write("%-9.3f" % hour_decimal[0])
outfile2.write("%-8.1f" % (Curr_sim_time/60.))
outfile2.write("%-8.1f" % Ins[0])
outfile2.write("%-8.1f" % Tamb[0])
outfile2.write("%-8.1f" % wind[0])

for value in Tg:
    outfile2.write("%-8.2f" % value)
for value in Tp:
    outfile2.write("%-8.2f" % value)    
for value in Tf:
    outfile2.write("%-8.2f" % value)    
for value in Tb:
    outfile2.write("%-8.2f" % value)    
outfile2.write("\n")
############################

Next_store_time += Store_dt                     # set next storage time
Next_store_time_details += Store_dt_details     # set next storage time in detailed output file

Tout_avg_stor_dt = 0.0           
Count_stor_dt = 1              # This counts how many timesteps have passed before storing results in the first output file, so that 
                               # the avg of outlet Tf during the store_dt can be calculated. Starts at 1 because the least is to write data every time step.

delta_Tg = ones(Nx+1)               # Matrix of iteration changes of Tg in each node, inititalized to a value above the tolerance desired (°C)
delta_Tp = ones(Nx+1)               # Matrix of iteration changes of Tp in each node, initialized to a value above the tolerance desired (°C)
delta_Tf = ones(Nx+1)               # Matrix of iteration changes of Tf in each node, inititalized to a value above the tolerance desired (°C)
delta_Tb = ones(Nx+1)               # Matrix of iteration changes of Tb in each node, initialized to a value above the tolerance desired (°C)


iterations_tot = 0    


# TIME LOOP STARTS HERE
while Curr_sim_time < Total_sim_time:
    iterations = 0
    
    h_dec = hour_decimal[int(Curr_sim_time/Data_time_step)]
    I = Ins[int(Curr_sim_time/Data_time_step)]
    Ta = Tamb[int(Curr_sim_time/Data_time_step)]
    v_wind = wind[int(Curr_sim_time/Data_time_step)]
    patm = Patm[int(Curr_sim_time/Data_time_step)]
    Xa = X[int(Curr_sim_time/Data_time_step)]
    Tdp = Tdewp[int(Curr_sim_time/Data_time_step)]
                
    mdot_f = mdot_a[int(Curr_sim_time/Data_time_step)]  # Air mass flow rate (kg/s)
    Vdot = mdot_f/rho_f                 # Airflow rate (m3/s)  
    v = Vdot/(W*H)                      # Average air speed in air channel (m/s)
    spec_Vdot = Vdot/Ac                 # Airflow rate per unit collector area (m3/s-m2)
    Re = v*Dh*rho_f/mu_f                # Reynolds number for flow in collector
    
    if Re <= 2550.0:                    # Calculates the Nusselt number based on flow conditions
        Nu = 5.385 + 0.148*Re*H/L
    elif Re <= 10000.0:
        Nu = 4.4e-4*Re**1.2 + 9.37*Re**0.471*H/L
    else:
        Nu = 0.03*Re**0.74 + 0.788*Re**0.74*H/L
        
    hc_fb = Nu*k_f/Dh               # Conv. heat transf. coeff. fluid-back plate (W/m2-K)
    hc_pf = Nu*k_f/Dh               # Conv. heat transf. coeff. absorber-fluid (W/m2-K)
    
    Ts = ((Ta + 273.15)*(0.711 + 0.0056*Tdp + 0.000073*Tdp**2 + 0.013*cos(15*h_dec*pi/180))**0.25) - 273.15 # Sky temperature (°C) (Duffie and Beckmann)
                                                                                                            # This equation makes use of the hour of the day (h_dec in the equation)
                                                                                                            # Here the h_dec is the hour of day in decimal, for example 6:30 am would be 6.5
                                                                                                            # and 2:15 pm would be 14.25            
    
    hc_gw = 5.7 + 3.8*v_wind        # Wind conv. heat trans. coeff. (W/m2-K)
    
    hr_pg = sigma*((Tp+273)**2 + (Tg+273)**2)*((Tp+273) + (Tg+273))/(1/e_p_up + 1/e_g - 1)        # Radiation heat transfer coefficients (W/m2-K)
    hr_pb = sigma*((Tp+273)**2 + (Tb+273)**2)*((Tp+273) + (Tb+273))/(1/e_p_down + 1/e_b - 1)
    hr_gs = sigma*e_g*((Tg+273)**2 + (Ts+273)**2)*((Tg+273) + (Ts+273))
    
    Ra = 9.81*(1/(Tf + 273))*abs(Tp - Tg)*(l**3)*rho_f/(mu_f*th_diff_f)         # Rayleigh number
    Nu_nc = 1 + 1.44*(1 - 1708*sin(1.8*tilt)**1.6/(Ra*cos(tilt)))*maximum(1 - 1708/(Ra*cos(tilt)), 0.0) + maximum((Ra*cos(tilt)/5830)**(1./3.) - 1, 0.0) # Nussel number for natural convection (Duffie and Beckman)
    
    hnc_pg = Nu_nc*k_f/l
    
    while (amax(delta_Tg) > TOL_delta_Tg) or (amax(delta_Tp) > TOL_delta_Tp) or (amax(delta_Tf) > TOL_delta_Tf) or (amax(delta_Tb) > TOL_delta_Tb):  # This is the iteration control, the tolerance value can be varied
        
        # LOOP OVER Tg NODES FROM i=0 TO Nx
        for j in range(0,Nx+1):
            
            jm1 = j-1 if j > 0 else j+1         # This is to be able to use the same equations for internal and boundary nodes
            jp1 = j+1 if j < Nx else j-1    
 
            Tg_iter = (Bg*Tg[j] + dt*(Ag*(Tg_new[jm1] + Tg_new[jp1]) + (hnc_pg[j] + hr_pg[j])*Tp_new[j] + alpha_g*I + hc_gw*Ta + hr_gs[j]*Ts))/(Bg + dt*(2*Ag + hnc_pg[j] + hc_gw + hr_pg[j] + hr_gs[j]))
            
            Tg_iter = relax*Tg_iter + (1-relax)*Tg_new[j]             # Applying underrelaxation 
            
            delta_Tg[j] = abs(Tg_iter - Tg_new[j])          # Recording the change in temperature from past iteration
                        
            Tg_new[j] = Tg_iter

        # LOOP OVER Tp NODES FROM i=0 TO Nx
        for j in range(0,Nx+1):
            
            jm1 = j-1 if j > 0 else j+1         # This is to be able to use the same equations for internal and boundary nodes
            jp1 = j+1 if j < Nx else j-1 

            Tp_iter = (Bp*Tp[j] + dt*((hnc_pg[j]+ hr_pg[j])*Tg_new[j] + Ap*(Tp_new[jm1] + Tp_new[jp1]) + hc_pf*Tf_new[j] + hr_pb[j]*Tb_new[j] + alpha_p*tau_g*I))/(Bp + dt*(2*Ap + hnc_pg[j] + hc_pf + hr_pg[j] + hr_pb[j]))
    
            Tp_iter = relax*Tp_iter + (1-relax)*Tp_new[j]             # Applying underrelaxation 
            
            delta_Tp[j] = abs(Tp_iter - Tp_new[j])              # Recording the change in temperature from past iteration
            
            Tp_new[j] = Tp_iter    
    
        # LOOP OVER AIR NODES FROM i=0 TO M                                
        Tf_iter = Ta
        Tf_new[0] = Tf_iter
        delta_Tf[0] = Tf_iter - Tf_new[0]
        for j in range(1,Nx+1):
            Tf_iter = (Tf[j] + dt*(Af*hc_pf*Tp_new[j] + (v/dx)*Tf_new[j-1] + Af*hc_fb*Tb_new[j]))/(1 + dt*(v/dx + Af*(hc_pf + hc_fb)))
            
            Tf_iter = relax*Tf_iter + (1-relax)*Tf_new[j]             # Applying underrelaxation 
            
            delta_Tf[j] = abs(Tf_iter - Tf_new[j])          # Recording the change in temperature from past iteration
            
            Tf_new[j] = Tf_iter
            
        # LOOP OVER Tb NODES FROM i=0 TO Nx
        for j in range(0,Nx+1):
            
            jm1 = j-1 if j > 0 else j+1         # This is to be able to use the same equations for internal and boundary nodes
            jp1 = j+1 if j < Nx else j-1 

            Tb_iter = (Bb*Tb[j] + dt*(Ub*Ta + hr_pb[j]*Tp_new[j] + hc_fb*Tf_new[j] + Ab*(Tb_new[jm1] +Tb_new[jp1])))/(Bb + dt*(2*Ab + hr_pb[j] + hc_fb + Ub))
        
            Tb_iter = relax*Tb_iter + (1-relax)*Tb_new[j]             # Applying underrelaxation 
            
            delta_Tb[j] = abs(Tb_iter - Tb_new[j])              # Recording the change in temperature from past iteration
            
            Tb_new[j] = Tb_iter
        
        iterations += 1
        iterations_tot += 1
    
    print ("max_delta_Tg", amax(delta_Tg))
    print ("max_delta_Tp", amax(delta_Tp))
    print ("max_delta_Tf", amax(delta_Tf))
    print ("max_delta_Tb", amax(delta_Tb))
     
    print ("t", Curr_sim_time)
    print ("iterations", iterations)    
    
    delta_Tg[:] = 1.0              # Reset the values of delta_H and delta_Ta in each node to a value above the tolerance desired
    delta_Tp[:] = 1.0
    delta_Tf[:] = 1.0
    delta_Tb[:] = 1.0
    
    Tg[:] = Tg_new                  # Replacing the new calculated temperatures into the vectors of current temperature to move into next timestep
    Tp[:] = Tp_new
    Tf[:] = Tf_new
    Tb[:] = Tb_new
    
    if Tf[Nx] > Ta:
        UsefulHeatTotal += (Tf[Nx]-Ta)*mdot_f*c_f*dt
        
    Total_I += I*Ac*dt
        
    Tout_avg_stor_dt =  (Tout_avg_stor_dt*(Count_stor_dt - 1) + Tf[Nx])/Count_stor_dt    # This is used to calculate the average outlet T during the store_dt.
                                                                                         # It calculates the mean of all Tout values between store_dts to give a value
                                                                                         # of mean Tout for the whole storage period. For example if i store data every minute
                                                                                         # and my ts is 10s, the code calculates a Tout every 10s, so here we calculate the mean
                                                                                         # of 6 Tout values which represent the mean Tout of that minute.                                                                                                                      

    Curr_sim_time += dt             # Advance time step
    
    
    # STORING DATA IN OUTPUT FILES
    
    index = int(Curr_sim_time/Data_time_step)
    if Curr_sim_time % Data_time_step == 0:
       index -= 1     
    
    if (Next_store_time - Curr_sim_time) < dt:
        outfile1.write("%-12s" % date[index])
        outfile1.write("%-9s" % hour_min[index])
        outfile1.write("%-9.3f" % hour_decimal[index])
        outfile1.write("%-8.1f" % (Curr_sim_time/60.))
        outfile1.write("%-8.1f" % Ins[index])
        outfile1.write("%-8.1f" % Tamb[index])
        outfile1.write("%-8.1f" % wind[index])
        outfile1.write("%-8.0f" % Patm[index])
        outfile1.write("%-8.4f" % X[index])
        outfile1.write("%-9.3f" % Tout_avg_stor_dt)
        outfile1.write("%-11.8f" % mdot_a[index])
        outfile1.write("\n")
        
        Tout_avg_stor_dt = 0
        Count_stor_dt = 0 
        Next_store_time += Store_dt
        
    if (Next_store_time_details - Curr_sim_time) < dt:
        
        outfile2.write("%-12s" % date[index])
        outfile2.write("%-9s" % hour_min[index])
        outfile2.write("%-9.3f" % hour_decimal[index])
        outfile2.write("%-8.1f" % (Curr_sim_time/60.))
        outfile2.write("%-8.2f" % Ins[index])
        outfile2.write("%-8.2f" % Tamb[index])
        outfile2.write("%-8.2f" % wind[index])
        
        for value in Tg:
            outfile2.write("%-8.2f" % value)
        for value in Tp:
            outfile2.write("%-8.2f" % value)    
        for value in Tf:
            outfile2.write("%-8.2f" % value)    
        for value in Tb:
            outfile2.write("%-8.2f" % value)
        if I > 0:
            outfile2.write("%-8.3f" % ((Tf[Nx]-Ta)*(mdot_f)*c_f/(I*Ac)))
        outfile2.write("\n")
       
        
        Next_store_time_details += Store_dt_details
    
    Count_stor_dt += 1
    

Col_eff = UsefulHeatTotal/Total_I

outfile2.write("Total I: %-15.2f" % Total_I)
outfile2.write("\n")
outfile2.write("Total Useful heat (J): %-15.2f" % UsefulHeatTotal)
outfile2.write("\n")
outfile2.write("Collector efficiency: %-15.3f" % Col_eff)
outfile2.write("\n")


outfile1.close()
outfile2.close() 
 
