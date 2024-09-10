# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 11:35:23 2022

@author: FranzR
"""
from numpy import log, zeros, linspace, ones, amax, exp, full, pi, mean, std
import sys

#########################
# HERE ARE TWO SUBROUTINES TO READ THE INPUT DATA FILE AND ASSIGN EACH DATA VARIABLE TO A LIST WHOSE LENGTH IS THE NUMBER OF DATA ROWS.
# AS IT IS, THE FIRST SUBROUTINE IS MADE TO READ THE OUPUT FILE FROM THE SOLAR AIR HEATER (SAH) PROGRAM WHICH HAS 6 ROWS WITH SETUP INFO, WHICH ARE SKIPPED:
# THE USER CAN MODIFY THE SUBROUTINE OR THE INPUT FILE TO MATCH IF THE DATA FILE FORMAT IS DIFFERENT
def extract_data_SAH_output(filename):
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
    Xa_in = []
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
        Xa_in.append(float(words[8]))
        Ta_in.append(float(words[9]))
        mdot_a.append(float(words[10]))
        
    infile.close()
        
    return date, hour_min, hour_dec, time_min, I, Tamb, wind, Patm, Xa_in, Ta_in, mdot_a

# THIS SECOND SUBROUTINE TO READ INPUT DATA IS MEANT TO USE THE OUTPUT DATA FROM THE LATENT HEAT STORAGE (LHS) AS INPUT. 
# IT DIFFERS FROM THE FIRST IN THAT THE LHS OUTPUT FILE HAS AN EXTRA COLUMN FOR THE LIQUID FRACTION WHICH IS NOT NEEDED HERE
# BUT IT HAS TO BE ACCOUNTED FOR
def extract_data_LHS_output(filename):
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
    Xa_in = []
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
        Xa_in.append(float(words[8]))
        Ta_in.append(float(words[9]))
        mdot_a.append(float(words[13]))
        
    infile.close()
        
    return date, hour_min, hour_dec, time_min, I, Tamb, wind, Patm, Xa_in, Ta_in, mdot_a
#####################

# FUNCTION TO CALCULATE THE SATURATION VAPOR PRESSURE FROM A GIVEN TEMPERATURE
def saturation_vapor_pressure(Tad):      # Saturation vapor pressure at some temperature, Pa (ASHRAE, 2017)
    Tabs = Tad + 273.15     
    log_Pvsat = -5.8002206e3/Tabs + 1.3914493e0 - 4.8640239e-2*Tabs + 4.1764768e-5*Tabs**2 - 1.4452093e-8*Tabs**3 + 6.5459673e0*log(Tabs)   
    Pvsat = exp(log_Pvsat)
    return Pvsat
###############################

# FUNCTION TO CALCULATE THE EQUILIBRIUM RELATIVE HUMIDITY BASED ON SORPTION ISOTHERM EQUATION
def equilibrium_rh(C1,C2,C3,Ta,Xp):
    rh_eq = (1./(((C1 + C2*Ta)/Xp)**C3 + 1.))
    return rh_eq

# FUNCTION TO CALCULATE THE EQUILIBRIUM MOISTURE CONTENT BASED ON SORPTION ISOTHERM EQUATION
def equilibrium_Xp(C1,C2,C3,Ta,rh):
    Xp_eq = (C1 + C2*Ta)*(rh/(1.-rh))**(1./C3)
    return Xp_eq

# FUNCTION TO CALCULATE THE K DRYING PARAMETER OF THE THIN LAYER DRYING EQUATION (HERE FOR THE PAGE MODEL)
def drying_rate_k_parameter(Ta,rh,v_a):
    k = 2.8e-3*exp(0.059*Ta)*(100*rh)**(-0.139)*v_a**0.025     # Page k parameter for wheat from Ramaj, 2021
    return k

# FUNCTION TO CALCULATE THE N DRYING PARAMETER OF THE THIN LAYER DRYING EQUATION (HERE FOR THE PAGE MODEL)
# IN THIS CASE N IS A CONSTANT BUT OFTEN IT IS NOT
def drying_rate_n_parameter():
    n = 0.784                                    # Page n parameter for wheat from Ramaj, 2021
    return n

# MAIN PROGRAM

Input_file = "LHS_Output_SAH_20m2_Slab3m_Stutt_7_8Aug2022 - SOFTX.txt"
Output_file = "FBD_Output_from_SAH20m2_LHS3m_Stutt_7_8Aug2022A.txt"   

# READING DATA FROM FILE
date, hour_min, hour_dec, t_min, Ins, Tamb, wind, P_atm, Xa_in, Ta_in, mdot_a = extract_data_LHS_output(Input_file)


# MAKE A RH_INLET LIST BASED ON THE INPUT DATA AIR CONDITIONS (THIS IS JUST TO BE ABLE TO WRITE THE INPUT RH TO THE OUTPUT FILE)
rh_in = []     
for i in range(len(Xa_in)):                 
                
    pv_Ta_out = Xa_in[i]*P_atm[i]/(0.621945 + Xa_in[i])
    pvsat_Ta_out = saturation_vapor_pressure(Ta_in[i])
    rh_in.append(pv_Ta_out/pvsat_Ta_out)
    
    
### USER INPUTS ###

hfg = 2253430.              # Latent heat of vaporization (asummed constant) (J/kg)
relax = 0.9                 # Underrelaxation factor (can be lower than 1 for stability if required)
Xp0_wetting = 0.0
Data_time_step = 60.    # This is the time resolution in seconds of the input data (s)

# INITIAL CONDITIONS
Mass = 1000.                # Mass of product at initial moisture content (kg)
Xp0 = 0.28205               # Initial product moisture content, db (kg/kg)
Tp0 = Tamb[0]               # Initial product temperature (assumed here to be the ambient temperature at t=0) (°C)
Ta0 = Tamb[0]               # Initial interstitial air temperature (assumed here to be the ambient temperature at t=0) (°C)

# PRODUCT PROPERTIES
C1 = 0.129                  # Parameters of Oswin isotherm equation for pionier WHEAT (paper Ramaj, 2021)
C2 = -6.46e-4               # This has to be changed by the user to parameters of the equation for his/her own product (there are various isotherm model equations)
C3 = 2.944                  # Care has to be taken not only with which isotherm equation is being used 
                            #  but also the units of the variables (rh as % or decimal, moisture as % or decimal and as wet basis or dry basis)

porosity = 0.4              # Porosity of the product bulk, from ASAE D241.4, -
rho_bulk = 774.4 - 703*(Xp0/(1+Xp0)) + 18510*(Xp0/(1+Xp0))**2 - 148960*(Xp0/(1+Xp0))**3 + 311600*(Xp0/(1+Xp0))**4    # Product bulk density depending on initial moisture (Nelson equation in ASAE D241.4 ) (kg/m3)
                                                                                                                     # This can be replaced by a single value    

ap = 1181.                  # Product specific volumetric area of product (Brooker, Bakker A. Table 2.6, 8.1 and D-5) (m2/m3)
c_dp = 1240.                # Product Specific heat of product in dry basis, from ASAE D243.4 (J/kg-K)                                               
Dp = 0.004                  # Product equivalent particle diameter for Reynolds number for heat trans. coeff. (from Brooker Tables 2.1 and 8.1) (m)

# AIR PROPERTIES
rho_a = 1.127               # Air density (at about 40 °C) (kg/m3)
c_a = 1007.                 # Specific heat of dry air (J/kg-K)
c_v = 1883.                 # Specific heat of water vapor (J/kg-K)
mu_a = 1.907e-5             # Air viscosity (Pa s)

c_w = 4187.                 # Specific heat of water (J/kg-K)

# DRYER DIMENSIONS
Side_length = 1.5            # Assumed that the dryer floor is of square shape (m)
Area_bin = Side_length**2    # Bin floor area (m2) (if rectangular bin)
#Dryer_diameter = 1.2                   # Dryer diameter (if round bin assumed) (m) 
#Area_bin = pi*(Dryer_diameter/2)**2    # Bin floor area (for round bin) (m2)

# SPACE AND TIME GRID SIZE
dt = 30.                          # Timestep (s)
dz = 0.008                       # Mesh size along dryer height (provisional) (m)

# TOLERANCE FOR THE DIFFERENCE BETWEEN CONSECUTIVE ITERATIONS
TOL_delta_Ta = 0.001            # Tolerance of nodal air temperature (°C)
TOL_delta_Xa = 0.00001          # Tolerance of nodal air humidity ratio (kg/kg)
TOL_delta_Tp = 0.001            # Tolerance of nodal product temperature (°C)
TOL_delta_Xp = 0.00001          # Tolerance of nodal product moisture content (kg/kg)

# CHOOSE HOW OFTEN TO STORE DATA
Store_dt_min = 1.    # Data storage interval in minutes (min)


# CALCULATIONS START HERE

Mp0 = Xp0/(1+Xp0)           # Initial product moisture content, wb (kg/kg)


rh0 = equilibrium_rh(C1,C2,C3,Ta0,Xp0)      # Relative humidity of interstitial air in equilibrium with product at beginning, -
pvsat0 = saturation_vapor_pressure(Ta0)     # Saturation vapor pressure at the initial interstitial temperature (Pa)
pv0 = rh0*pvsat0                            # Actual initial vapor pressure of the interstitial air (Pa)
Xa0 = 0.621945*pv0/(P_atm[0] - pv0)         # Initial humidity ratio of interstitial air (kg/kg)

rho_dp_bulk = rho_bulk*(1-Mp0)              # Product bulk dry matter density (kg/m3)
                                            # Since the volume of the grain bed is fixed in the model the dry matter bulk density is constant
Bulk_volume = Mass/rho_bulk                 #Volume of bulk (m3) 
                                    

Height = Bulk_volume/Area_bin               # Height of the grain bed (m) 
                                                                                
# CREATING THE NECESSARY VARIABLES, VECTORS AND ARRAYS FOR THE SCHEME
M = int(Height//dz)           # Number of divisions along dryer height
numb_stor_div = 8             # This is the number of divisions along the bed height in which data should be store, first position being bottom and last being top of the bed 
    
while M%(numb_stor_div) != 0:       # This is to find a number of divisions for the bed that is multiple of numb_stor_div
    M +=1                            
    
dz = Height/M                             # Mesh size (actual) (m)
z_position = linspace(0, Height, M+1)     # Vector of node coordinates along the dryer height (m)        

write_pos_index = []                                # Here a list of indexes of the z_position vector is created where we want to store results
for i in range(0, numb_stor_div+1):
    write_pos_index.append(int(i/(numb_stor_div)*M))
       
Ta = full(M+1, Ta0)                 # Vector of interstitial air temperature along dryer height (°C)
Ta_new = full(M+1, Ta0)             # Vector of interstitial air temperature along dryer height at new time level (°C)

Xa = full(M+1, Xa0)                 # Vector of interstitial air humidity along dryer height (°C)
Xa_new = full(M+1, Xa0)             # Vector of interstitial air humidity along dryer height at new time level (°C)

Tp = full(M+1, Tp0)                 # Vector of product temperature along dryer height (°C)
Tp_new = full(M+1, Tp0)             # Vector of product temperature along dryer height in the next time level (°C)

Xp = full(M+1, Xp0)                 # Vector of product moisture along dryer height (°C)
Xp_new = full(M+1, Xp0)             # Vector of product moisture along dryer height in the next time level (°C)


Store_dt = Store_dt_min*60  # Data storage interval in seconds (s)
Total_sim_t = len(t_min)*Data_time_step      # The total simulation time (s)
     
delta_Ta = ones(M+1)               # Matrix of iteration changes of Ta in each node, inititalized to a value above the tolerance desired
delta_Xa = ones(M+1)                 # Matrix of iteration changes of Xa in each node, initialized to a value above the tolerance desired
delta_Tp = ones(M+1)               # Matrix of iteration changes of Tp in each node, inititalized to a value above the tolerance desired
delta_Xp = ones(M+1)                 # Matrix of iteration changes of Xp in each node, initialized to a value above the tolerance desired

iterations_tot = 0    
Curr_sim_time = 0.0                       # The current simulation time during the calculations (s)
Next_store_time = 0.

# WRITE SETUP INFO AND OUTPUT DATA TO FILE
outfile = open(Output_file, "w")
outfile.write( "Mass:%-8.1f Volume:%-8.4f Xp0:%-6.1f Mp0:%-6.2f Tp0:%-5.2f porosity:%-5.3f" % (Mass, Bulk_volume, Xp0, Mp0, Tp0, porosity))
outfile.write("\n")
outfile.write( "rho_bulk:%-6.2f rho_dp_bulk:%-6.2f ap:%-8.2f c_dp:%-8.2f Dp:%-7.4f" % (rho_bulk, rho_dp_bulk, ap, c_dp, Dp))
outfile.write("\n")
outfile.write( "Side_length:%-5.2f Bin_area:%-6.3f Bulk_height:%-5.3f Ta0:%-5.2f Xa_0:%-6.4f" % (Side_length, Area_bin, Height, Ta0, Xa0))
outfile.write("\n")
outfile.write( "dt:%-6.1f dz:%-6.4f" % (dt, dz))
outfile.write("\n")

outfile.write("%-12s%-9s%-9s%-8s%-8s%-8s%-8s%-8s%-8s%-8s%-8s"  % ("Date", "Hour_min", "Hour_dec", "t_min", "Ins", "Tamb", "Wind", "Patm", "Xa_in", "Ta_in", "rh_in" ))

for ind in write_pos_index:
    outfile.write("%-8.4f" % z_position[ind])  
outfile.write("%-8s" % ("Xp_avg"))
outfile.write("%-8s" % ("Xp_std"))

for ind in write_pos_index:
    outfile.write("%-8.4f" % z_position[ind])    
outfile.write("%-8s"  % ("Tp_avg"))

for ind in write_pos_index:
    outfile.write("%-8.4f" % z_position[ind])    
outfile.write("%-8s" % ("Xa_avg"))

for ind in write_pos_index:
    outfile.write("%-8.4f" % z_position[ind])    
outfile.write("%-8s" % ("Ta_avg"))

outfile.write("\n")

outfile.write("%-12s" % date[0])
outfile.write("%-9s" % hour_min[0])
outfile.write("%-9s" % hour_dec[0])
outfile.write("%-8.1f" % (Curr_sim_time/60.))
outfile.write("%-8.1f" % Ins[0])
outfile.write("%-8.1f" % Tamb[0])
outfile.write("%-8.1f" % wind[0])
outfile.write("%-8.0f" % P_atm[0])
outfile.write("%-8.4f" % Xa_in[0])
outfile.write("%-8.2f" % Ta_in[0])
outfile.write("%-8.3f" % rh_in[0])

for ind in write_pos_index:
    outfile.write("%-8.4f" % Xp[ind])

outfile.write("%-8.3f" % mean(Xp))
outfile.write("%-8.5f" % std(Xp))

for ind in write_pos_index:
    outfile.write("%-8.2f" % Tp[ind])

outfile.write("%-8.2f" % mean(Tp))

for ind in write_pos_index:
    outfile.write("%-8.4f" % Xa[ind])

outfile.write("%-8.4f" % mean(Xa))

for ind in write_pos_index:
    outfile.write("%-8.2f" % Ta[ind])

outfile.write("%-8.2f" % mean(Ta))

outfile.write("\n")
#########################

Next_store_time += Store_dt


# TIME LOOP STARTS HERE
while Curr_sim_time < Total_sim_t:

    iterations = 0
    Tair_in = Ta_in[int(Curr_sim_time/Data_time_step)]      # Get current inlet air temperature (°C)
    Xair_in = Xa_in[int(Curr_sim_time/Data_time_step)]      # Get humidity ratio of inlet air (kg/kg)
    patm = P_atm[int(Curr_sim_time/Data_time_step)]         # Get atmospheric pressure (Pa)
    mdot_air = mdot_a[int(Curr_sim_time/Data_time_step)]    # Get current air mass flow rate (kg/s)
    V_a = mdot_air/rho_a                                    # Airflow rate (m3/s)
    v_a = V_a/Area_bin                                      # Superficial velocity of drying air (m/s)
    Re = v_a*rho_a*Dp/mu_a                                  # Reynolds number
    
    h = 0.2755*(c_a/1000.)*(v_a*3600)*rho_a*Re**(-0.34)     # Heat transfer coefficient according to Brooker, B-A book and paper de B-A 1967
    A = dt*v_a/(porosity*dz)                                # Common factor in eqs for Ta, Xa
    
    while (amax(delta_Ta) > TOL_delta_Ta) or (amax(delta_Xa) > TOL_delta_Xa) or (amax(delta_Tp) > TOL_delta_Tp) or (amax(delta_Xp) > TOL_delta_Xp):  # This is the iteration control, the tolerance value can be varied
        
        # LOOP OVER Ta NODES FROM i=0 TO M                         
        Ta_iter = Tair_in
        Ta_new[0] = Ta_iter
        delta_Ta[0] = Ta_iter - Ta_new[0]
        for j in range (1,M+1):
            
            rho_da = 1./(0.287042*(Ta_new[j] + 273.15)*(1 + 1.607858*Xa_new[j])/(0.001*patm)) # Calculate dry air density for current node. First trying with Ta[j] and Xa[j] instead of new values.
            if rho_da != rho_da:
                sys.exit() 
             
            B = porosity*rho_da*(c_a + c_v*Xa_new[j]) # Repeated denominator of eq for Ta
            if B != B:
                sys.exit()
                            
            Ta_iter = (Ta[j] + A*Ta_new[j-1] + dt*h*ap*Tp_new[j]/B)/(1. + A + dt*h*ap/B)   
            if Ta_iter != Ta_iter:
                sys.exit()                                                                                         
                                                                                                    
            Ta_iter = relax*Ta_iter + (1-relax)*Ta_new[j]             # Applying underrelaxation 
            delta_Ta[j] = abs(Ta_iter - Ta_new[j])                    # Recording the change in temperature from past iteration
            Ta_new[j] = Ta_iter            
            
        # LOOP OVER Xa NODES FROM i=0 TO M
        Xa_iter = Xair_in
        Xa_new[0] = Xa_iter
        delta_Xa[0] = Xa_iter - Xa_new[0]
        
        for j in range(1,M+1):
                    
            rho_da = 1/(0.287042*(Ta_new[j] + 273.15)*(1 + 1.607858*Xa_new[j])/(0.001*patm)) # Calculate dry air density for current node
            
            if rho_da != rho_da:
                sys.exit()
                
            Xa_iter = (Xa[j] + A*Xa_new[j-1] - rho_dp_bulk/(porosity*rho_da)*(Xp_new[j] - Xp[j]))/(1 + A)   # Here also using Xp upwind values since air humidity will be
            if Xa_iter != Xa_iter:                                                                          # more influenced by the product moisture content upwind of the air
                sys.exit()                                                                                  # node than by the downwind product node                  
                                                                                                                 
            Xa_iter = relax*Xa_iter + (1-relax)*Xa_new[j]             # Applying underrelaxation 
            
            delta_Xa[j] = abs(Xa_iter - Xa_new[j])                  # Recording the change in humidity ratio from past iteration
            
            Xa_new[j] = Xa_iter
            
            
        # LOOP OVER Tp NODES FROM i=0 TO M 
        j = 0
        rho_da = 1/(0.287042*(Ta_new[j] + 273.15)*(1 + 1.607858*Xa_new[j])/(0.001*patm)) # Calculate dry air density for current node
        if rho_da != rho_da:
            sys.exit()
                    
        C = (1 + (h*ap*dt*dz - rho_da*v_a*c_v*dt*(Xa_new[j+1] - Xa_new[j]))/(rho_dp_bulk*dz*(c_dp + c_w*Xp_new[j])))
        if C != C:
            sys.exit()
            
        Tp_iter = (Tp[j] + (dt/(rho_dp_bulk*(c_dp + c_w*Xp_new[j])))*(h*ap*Ta_new[j] - (hfg + c_v*Ta_new[j])*(rho_da*v_a*(Xa_new[j+1] - Xa_new[j])/dz)))/C
        if Tp_iter != Tp_iter:
            sys.exit()
        
        Tp_iter = relax*Tp_iter + (1-relax)*Tp_new[j]             # Applying underrelaxation  
        
        delta_Tp[j] = abs(Tp_iter - Tp_new[j])                  # Recording the change in temperature from past iteration
        
        Tp_new[j] = Tp_iter
            
        for j in range(1,M):
            
            rho_da = 1/(0.287042*(Ta_new[j] + 273.15)*(1 + 1.607858*Xa_new[j])/(0.001*patm)) # Calculate dry air density for current node
            if rho_da != rho_da:
                sys.exit()
            
            C = (1 + (2*h*ap*dt*dz - rho_da*v_a*c_v*dt*(Xa_new[j+1] - Xa_new[j-1]))/(2*rho_dp_bulk*dz*(c_dp + c_w*Xp_new[j])))
            if C != C:
                sys.exit()
            
            Tp_iter = (Tp[j] + (dt/(rho_dp_bulk*(c_dp + c_w*Xp_new[j])))*(h*ap*Ta_new[j] - (hfg + c_v*Ta_new[j])*(rho_da*v_a*(Xa_new[j+1] - Xa_new[j-1])/(2*dz))))/C
            if Tp_iter != Tp_iter:
                sys.exit()
                
            Tp_iter = relax*Tp_iter + (1-relax)*Tp_new[j]             # Applying underrelaxation    
        
            delta_Tp[j] = abs(Tp_iter - Tp_new[j])                  # Recording the change in temperature from past iteration
        
            Tp_new[j] = Tp_iter
            
        j = M
        rho_da = 1/(0.287042*(Ta_new[j] + 273.15)*(1 + 1.607858*Xa_new[j])/(0.001*patm)) # Calculate dry air density for current node 
        if rho_da != rho_da:
            sys.exit()
                    
        C = (1 + (h*ap*dt*dz - rho_da*v_a*c_v*dt*(Xa_new[j] - Xa_new[j-1]))/(rho_dp_bulk*dz*(c_dp + c_w*Xp_new[j])))
        if C != C:
            sys.exit()
            
        Tp_iter = (Tp[j] + (dt/(rho_dp_bulk*(c_dp + c_w*Xp_new[j])))*(h*ap*Ta_new[j] - (hfg + c_v*Ta_new[j])*(rho_da*v_a*(Xa_new[j] - Xa_new[j-1])/dz)))/C
        if Tp_iter != Tp_iter:
            sys.exit()
        
        Tp_iter = relax*Tp_iter + (1-relax)*Tp_new[j]             # Applying underrelaxation
        
        delta_Tp[j] = abs(Tp_iter - Tp_new[j])                  # Recording the change in temperature from past iteration
        
        Tp_new[j] = Tp_iter  
            
        # LOOP OVER PRODUCT MOISTURE CONTENT FROM i=0 TO M-1
        # HERE IT IS MADE USE OF A THIN LAYER DRYING EQUATION WHICH IS PRODUCT SPECIFIC
        # AS WITH SORPTION ISOTHERMS, THERE ARE NUMEROUS MODELS. ONE USUALLY ADEQUATE FOR AGRIC. PRODUCTS IS THE PAGE MODEL; WHICH IS USED BELOW
        # IF THE SIMPLRE LEWIS MODEL IS TO BE USED, THE PARAMETER N SHOULD BE SET TO 1.
        # AS WITH THE ISOTHERM MODEL, CARE MUST BE TAKEN TO ENSURE THAT THE UNITS ARE TRANSFORMED TO THE UNITS OF THE PROGRAM IF REQUIRED.
        # MOISTURE CONTENT SHOULD BE IN DRY BASIS, TIME IN SECONDS
        # AS IT IS, THE PROGRAM USES THE SAME DRYING EQUATION FOR REWETTING. IF A REWETTING EQUATION IS APPROPRIATE THE USER CAN CAREFULLY MODIFY THE CODE
            
        for j in range(0,M+1):
            
            pvsat_Ta = saturation_vapor_pressure(Ta_new[j])  # Saturation vapor pressure at the air temperature in the current node (Pa)
            if pvsat_Ta != pvsat_Ta:
                sys.exit()
                
            pv_Ta = Xa_new[j]*patm/(0.621945 + Xa_new[j])          # Actual vapor pressure of air in current node (Pa)
            if pv_Ta != pv_Ta:
                sys.exit()
            
            rh = pv_Ta/pvsat_Ta 
            if rh != rh:
                sys.exit()
                              
            if rh >= 0.99:
                rh = 0.99
            Xp_eq = equilibrium_Xp(C1,C2,C3,Ta_new[j],rh)
            if Xp_eq != Xp_eq:
                sys.exit()
            
                
            k = drying_rate_k_parameter(Ta_new[j],rh,v_a)
            n = drying_rate_n_parameter()
            
                        
            Xp_current = Xp_new[j]
            Xp_init = Xp0
                        
            if Xp_current > Xp_eq:          # In this case drying should happen in that layer. We assume that Xp_current is never exactly the X_eq
                if Xp_eq == Xp_init:
                    Xp_eq = Xp_init - 0.0001
                if Xp_current > Xp_init:
                    Xp_init = Xp_current
                MR = (Xp_current-Xp_eq)/(Xp_init - Xp_eq)
            
                teq = (-log(MR)/k)**(1./n)*60.
                       
            else:
                MR = (Xp_current-Xp_eq)/(Xp0_wetting - Xp_eq)    # In case Xp_current is less than Xp_eq, wetting should happen, and we use a special initial moisture content
                if MR != MR:
                    sys.exit()
                    
                teq = (-log(MR)/k)**(1./n)*60.          # Multiplied by 60 because the drying equation was fitted with time in min but this model works in seconds
            
            Xp_iter = (Xp[j] + (dt/60.)*n*k*((teq+dt)/60.)**(n-1)*Xp_eq)/(1 + (dt/60.)*n*k*((teq+dt)/60.)**(n-1))
            if Xp_iter != Xp_iter:
                sys.exit()
            Xp_iter = relax*Xp_iter + (1-relax)*Xp_new[j]             # Applying underrelaxation  
            
            delta_Xp[j] = abs(Xp_iter - Xp_new[j])                   # Recording the change in moisture content from past iteration
            
            
            Xp_new[j] = Xp_iter
            
        iterations += 1
        iterations_tot += 1
    
    
    print ("max_delta_Ta", amax(delta_Ta))
    print ("max_delta_Xa", amax(delta_Xa))
    print ("max_delta_Tp", amax(delta_Tp))
    print ("max_delta_Xp", amax(delta_Xp))
    print ("t", Curr_sim_time)
    print ("iterations", iterations)       
    
    delta_Ta[:] = 1.0              # Reset the values of the deltas in each node to a value above the tolerance desired
    delta_Xa[:] = 1.0
    delta_Tp[:] = 1.0
    delta_Xp[:] = 1.0
    
    Curr_sim_time += dt             # Advance time step
    
    
    Ta[:] = Ta_new                  # Replacing the new calculated values into the vectors of current temperature to move into next timestep
    Xa[:] = Xa_new
    Tp[:] = Tp_new
    Xp[:] = Xp_new
    
    index = int(Curr_sim_time/Data_time_step)
    if Curr_sim_time % Data_time_step == 0:
       index -= 1
    
    # STORING DATA IN OUTPUT FILE  
    
    if (Next_store_time - Curr_sim_time) < dt:
        
        outfile.write("%-12s" % date[index])
        outfile.write("%-9s" % hour_min[index])
        outfile.write("%-9s" % hour_dec[index])
        outfile.write("%-8.1f" % (Curr_sim_time/60.))
        outfile.write("%-8.1f" % Ins[index])
        outfile.write("%-8.1f" % Tamb[index])
        outfile.write("%-8.1f" % wind[index])
        outfile.write("%-8.0f" % P_atm[index])
        outfile.write("%-8.4f" % Xa_in[index])
        outfile.write("%-8.2f" % Ta_in[index])
        outfile.write("%-8.3f" % rh_in[index])
        
        for ind in write_pos_index:
            outfile.write("%-8.4f" % Xp[ind])
        
        outfile.write("%-8.3f" % mean(Xp))
        outfile.write("%-8.5f" % std(Xp))
        
        for ind in write_pos_index:
            outfile.write("%-8.2f" % Tp[ind])

        outfile.write("%-8.2f" % mean(Tp))

        for ind in write_pos_index:
            outfile.write("%-8.4f" % Xa[ind])

        outfile.write("%-8.4f" % mean(Xa))

        for ind in write_pos_index:
            outfile.write("%-8.2f" % Ta[ind])

        outfile.write("%-8.2f" % mean(Ta))
        
        outfile.write("\n")

                
        Next_store_time += Store_dt
       
outfile.close()

        