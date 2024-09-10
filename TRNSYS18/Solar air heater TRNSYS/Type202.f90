  SUBROUTINE TYPE202

! USE STATEMENTS 
Use TrnsysFunctions 
!-----------------------------------------------------------------------

!export this subroutine for its use in external DLLs. 
!DEC$ATTRIBUTES DLLEXPORT :: TYPE202
    
implicit none
    
! LOCAL VARIABLE DECLARATIONS
    
!   From Trnsys
double precision :: Time, dt_hours ! the current simulation time and the timestep (h)
 

!   What are considered TRNSYS parameters, meaning user inputs that dont change        
double precision :: relax, L, W, H, gap_gp, tilt_deg                                                ! 6
double precision :: thick_g, thick_p, thick_b, k_g, k_p, k_b, c_g, c_p, c_b, rho_g, rho_p, rho_b    ! 12
double precision :: thick_i, k_i, alpha_g, alpha_p, alpha_b, tau_g, e_g, e_p_up, e_p_down, e_b      ! 10
double precision :: rho_f, c_f, mu_f, k_f                                                           ! 4
double precision :: dx_tentative,T0                                                                 ! 2
double precision :: TOL_delta_Tg, TOL_delta_Tp, TOL_delta_Tf, TOL_delta_Tb                          ! 4
    
integer, parameter :: Nx_max = 100                            ! 1, Maximum number of x nodes along collector length
double precision :: sigma, pi

!   What are considered TRNSYS inputs, meaning time dependent user inputs
double precision :: h_dec, Irr, Ta, rh, v_wind, p_atm, mdot_f, Tf_in  

!   Internal variables necessary which only need to be calculated once
double precision :: Dh, Ub, th_diff_f, Ac, tilt_rad
double precision :: Ag, Bg, Ap, Bp, Af, Ab, Bb
double precision :: dx, dt
integer :: Nx, Static_storage_size
double precision :: x_position(0:Nx_max)

!   Internal variables necessary which may need to be calculated every iteration or timestep
double precision :: Xa, Tdp, Ts, Vdot
double precision :: v, spec_Vdot, Re
double precision :: Nu, hc_fb, hc_pf, hc_gw
double precision :: Tg_iter, Tp_iter, Tf_iter, Tb_iter, Variable_update, eff, Usef_Q_rate, iterations
integer :: i, j, jm1, jp1
double precision :: Tabs, Pvsat, Pv   ! Variables to calculate humidity ratio and Tdp
double precision :: hr_pg(0:Nx_max), hr_gs(0:Nx_max), hr_pb(0:Nx_max), Ra(0:Nx_max), Nu_nc(0:Nx_max), hnc_pg(0:Nx_max)    
double precision :: Tg(0:Nx_max), Tg_new(0:Nx_max), Tp(0:Nx_max), Tp_new(0:Nx_max), Tf(0:Nx_max), Tf_new(0:Nx_max), Tb(0:Nx_max), Tb_new(0:Nx_max)
double precision :: delta_Tg(0:Nx_max), delta_Tp(0:Nx_max), delta_Tf(0:Nx_max), delta_Tb(0:Nx_max) 
integer :: CurrentUnit, CurrentType                         ! used in Message calls if there is an error (from TRNSYS forum)    
    

Data sigma/5.6697e-8/, pi/3.1415927/
!---------------------------------------------------------------------- 
!Get the Global TNRSYS Simulation Variables 
Time = getSimulationTime() 
dt_hours = getSimulationTimeStep()
CurrentUnit = getCurrentUnit()
CurrentType = getCurrentType()

!------------------------------------------------------------------------------
! VERSION SIGNING 
If (getIsVersionSigningTime()) Then 
    Call SetTypeVersion(17) 
    Return 
Endif 

!------------------------------------------------------------------------------
!DO ALL OF THE "LAST CALL" MANIPULATIONS HERE (AS IT IS IT DOESNT DO ANYTHING. MAYBE CHANGE IT TO CALCULATE SOME AGGREGATE RESULT) 
If (getIsLastCallofSimulation()) Then
    If (getNumberOfErrors() > 0) Then
        print *, "Terminated due to error"
    Endif
    Return
Endif

!---------------------------------------------------------------------- 
!   PERFORM ANY "END OF TIME STEP" MANIPULATIONS THAT MAY BE REQUIRED (AQUI UPDATE STORAGE VARIABLES) 
If (getIsEndOfTimestep()) Then
    L               = getParameterValue(2)   ! Collector length (m)
    dx_tentative    = getParameterValue(33)   ! tentative (might differ from definitive) grid size along collector length (m)       
    
    Nx = nint(L/dx_tentative)               ! Number of divisions along collector length (x direction), -    
    
    Static_storage_size = 8*(Nx+1)           ! Total static storage spots needed for T's and T_new's
    
    !   UPDATE STORED VALUES FROM THE T_new's to the T's
    do j = Static_storage_size/2+1, Static_storage_size 
        Variable_update = getStaticArrayValue(j)
        call setStaticArrayValue(j-Static_storage_size/2, Variable_update)
    end do
    Return 
Endif

!---------------------------------------------------------------------- 
!   DO ALL OF THE "INITIALIZATION CALL  MANIPULATIONS HERE 
If (getIsFirstCallofSimulation()) Then
    L               = getParameterValue(2)   ! Collector length (m)
    dx_tentative    = getParameterValue(33)   ! tentative (might differ from definitive) grid size along collector length (m)       
    Nx = nint(L/dx_tentative)               ! Number of divisions along collector length (x direction), -    
    
    Static_storage_size = 8*(Nx+1)           ! Total static storage spots needed for T's and T_new's


!   Tell the TRNSYS Engine How This Type Works 
    Call SetNumberofParameters(38)            
    Call SetNumberofInputs(8)            
    Call SetNumberofDerivatives(0)            
    Call SetNumberofOutputs(6)            
    Call SetIterationMode(1)     ! 1 indicates the type should be caller every iteration       
    Call SetNumberStoredVariables(Static_storage_size,0) 
    
    !  Set the Correct Input and Output Variable Types 
    Call SetInputUnits(1,'DM1')     ! Solar time (hours decimal)
    Call SetInputUnits(2,'IR2')     ! Solar radiation (W/m2)
    Call SetInputUnits(3,'TE1')     ! Ambient temperature (degC)
    Call SetInputUnits(4,'DM1')     ! Relative humidity (-)
    Call SetInputUnits(5,'VE1')     ! Wind speed (m/s)
    Call SetInputUnits(6,'PR3')     ! Atmospheric pressure (Pa)
    Call SetInputUnits(7,'MF2')     ! Mass flow rate (kg/s)
    Call SetInputUnits(8,'TE1')     ! Inlet temperature (degC)
  
    
    Call SetOutputUnits(1,'TE1')    ! Outlet air temperature (degC)    
    Call SetOutputUnits(2,'MF2')    ! Mass flow rate (kg/s)
    Call SetOutputUnits(3,'DM1')    ! Collector efficiency (-)
    Call SetOutputUnits(4,'PW2')    ! Useful energy rate (W)
    Call SetOutputUnits(5,'DM1')    ! Number of iterations 
    Call SetOutputUnits(6,'DM1')    ! Outlet air humidity ratio from ambient conditions (-)
    Return 
EndIf 
    
    
 !------------------------------------------------------------------------------ 
    !DO ALL OF THE "START TIME" MANIPULATIONS HERE - There Are No Iterations at the Intial Time 
If (getIsStartTime()) Then
    !  Read in the Values of the Parameters from the Input File 
    relax       = getParameterValue(1) ! Underrelaxation factor (between 0 and 1)
    L           = getParameterValue(2) ! Collector length (m) 
    W           = getParameterValue(3) ! Collector width (m)
    H           = getParameterValue(4) ! Air channel height (m)
    gap_gp      = getParameterValue(5) ! Gap bewtween cover and absorber (m)
    tilt_deg    = getParameterValue(6) ! Collector tilt from horizontal (deg)
    thick_g     = getParameterValue(7) ! Cover thickness (m)
    thick_p     = getParameterValue(8) ! Absorber thickness (m)
    thick_b     = getParameterValue(9) ! Bottom plate thickness (m)
    k_g         = getParameterValue(10) ! Cover conductivity (W/m-K)
    k_p         = getParameterValue(11) ! Absorber conductivity (W/m-K)
    k_b         = getParameterValue(12) ! Bottom plate conductivity (W/m-K)
    c_g         = getParameterValue(13) ! Cover specific heat (J/kg-K)
    c_p         = getParameterValue(14) ! Absorber specific heat (J/kg-K)
    c_b         = getParameterValue(15) ! Bottom plate specific heat (J/kg-K)
    rho_g       = getParameterValue(16) ! Cover density (kg/m3)
    rho_p       = getParameterValue(17) ! Absorber density (kg/m3)
    rho_b       = getParameterValue(18) ! Bottom plate density (kg/m3)
    thick_i     = getParameterValue(19) ! Insulation thickness (m)
    k_i         = getParameterValue(20) ! Insulation conductivity (W/m-K)
    alpha_g     = getParameterValue(21) ! Cover absorptance (-) 
    alpha_p     = getParameterValue(22) ! Absorber absorptance (-) 
    alpha_b     = getParameterValue(23) ! Bottom plate absorptance (-)
    tau_g       = getParameterValue(24) ! Cover transmittance (-)
    e_g         = getParameterValue(25) ! Cover emittance (-)
    e_p_up      = getParameterValue(26) ! Absorber emittance upwards (-)
    e_p_down    = getParameterValue(27) ! Absorber emittance downwards (-)
    e_b         = getParameterValue(28) ! Bottom plate emittance (-) 
    rho_f       = getParameterValue(29) ! Air density (considered constant) (kg/m3)
    c_f         = getParameterValue(30) ! Air specific heat (considered constant) (J/kg-K)
    mu_f        = getParameterValue(31) ! Air dynamic viscosity (considered constant) (Pa-s)
    k_f         = getParameterValue(32) ! Air thermal conductivity (considered constant) (W/m-K) 
    dx_tentative    = getParameterValue(33) ! tentative (might differ from definitive) grid size along collecto length (m) 
    TOL_delta_Tg    = getParameterValue(34) ! tolerance for error between iterations for cover temperature (degC)
    TOL_delta_Tp    = getParameterValue(35) ! tolerance for error between iterations for absorber temperature (degC)
    TOL_delta_Tf    = getParameterValue(36) ! tolerance for error between iterations for fluid temperature (degC)
    TOL_delta_Tb    = getParameterValue(37) ! tolerance for error between iterations for bottom plate temperature (degC)
    T0              = getParameterValue(38) ! Initial collector temperature (degC)

   !   Check the parameters for problems and RETURN if any are found 
    IF((relax .GT. 1.).OR.(relax .LE. 0.)) CALL foundBadParameter(1,'Fatal','The underrelax factor must be above 0 and up to 1')
    IF(L .LE. 0.)         CALL foundBadParameter(2,'Fatal','The collector length must be positive') 
    IF(W .LE. 0.)         CALL foundBadParameter(3,'Fatal','The collector width must be positive')
    IF(H .LE. 0.)         CALL foundBadParameter(4,'Fatal','The air channel height must be positive') 
    IF(gap_gp .LE. 0.)   CALL foundBadParameter(5,'Fatal','The gap between cover and absorber must be positive')
    IF((tilt_deg .LT. 0.).OR.(tilt_deg .GT. 90))      CALL foundBadParameter(6,'Fatal','The collector tilt angle must be between 0 and 90 deg')
    IF(thick_g .LE. 0.)   CALL foundBadParameter(7,'Fatal','The cover thicknesst must be positive')
    IF(thick_p .LE. 0.)   CALL foundBadParameter(8,'Fatal','The absorber thickness must be positive') 
    IF(thick_b .LE. 0.)   CALL foundBadParameter(9,'Fatal','The bottom plate must be positive')
    IF(k_g .LE. 0.)       CALL foundBadParameter(10,'Fatal','The cover conductivity must be positive')
    IF(k_p .LE. 0.)       CALL foundBadParameter(11,'Fatal','The absorber conducticity must be positive')
    IF(k_b .LE. 0.)       CALL foundBadParameter(12,'Fatal','The bottom plate conductivity must be positive')        
    IF(c_g .LE. 0.)       CALL foundBadParameter(13,'Fatal','The cover specific heat must be positive')
    IF(c_p .LE. 0.)       CALL foundBadParameter(14,'Fatal','The absorber specific heat must be positive')
    IF(c_b .LE. 0.)       CALL foundBadParameter(15,'Fatal','The bottom plate specific heat must be positive')
    IF(rho_g .LE. 0.)     CALL foundBadParameter(16,'Fatal','The cover density must be positive')
    IF(rho_p .LE. 0.)     CALL foundBadParameter(17,'Fatal','The absorber density must be positive')
    IF(rho_b .LE. 0.)     CALL foundBadParameter(18,'Fatal','The bottom plate density must be positive')
    IF(thick_i .LE. 0.)   CALL foundBadParameter(19,'Fatal','The back insulation thickness must be positive')
    IF(k_i .LE. 0)        CALL foundBadParameter(20,'Fatal','The back insulation conductivity must be positive')
    IF((alpha_g .LE. 0.).OR.(alpha_g .GT. 1.0))     CALL foundBadParameter(21,'Fatal','The cover absorptance must be above 0 and up to 1')
    IF((alpha_p .LE. 0.).OR.(alpha_p .GT. 1.0))     CALL foundBadParameter(22,'Fatal','The absorber absorptance must be above 0 and up to 1')
    IF((alpha_b .LE. 0.).OR.(alpha_b .GT. 1.0))     CALL foundBadParameter(23,'Fatal','The bottom plate absorptance must be above 0 and up to 1')
    IF((tau_g .LE. 0.).OR.(tau_g .GT. 1.0))         CALL foundBadParameter(24,'Fatal','The cover transmittance must be above 0 and up to 1')
    IF((e_g .LT. 0.).OR.(e_g .GT. 1.0))             CALL foundBadParameter(25,'Fatal','The cover emittance must be between 0 and 1')
    IF((e_p_up .LT. 0.).OR.(e_p_up .GT. 1.0))       CALL foundBadParameter(26,'Fatal','The absorber emittance upwards must be between 0 and 1')
    IF((e_p_down .LT. 0.).OR.(e_p_down .GT. 1.0))   CALL foundBadParameter(27,'Fatal','The absorber emittance downwards must be between 0 and 1')
    IF((e_b .LT. 0.).OR.(e_b .GT. 1.0))             CALL foundBadParameter(28,'Fatal','The bottom plate emittance must be between 0 and 1')
    IF(rho_f .LE. 0.)     CALL foundBadParameter(29,'Fatal','The fluid density must be positive')
    IF(c_f .LE. 0.)       CALL foundBadParameter(30,'Fatal','The fluid specific heat must be positive')
    IF(mu_f .LE. 0.)      CALL foundBadParameter(31,'Fatal','The fluid dynamic viscosity must be positive')
    IF(k_f .LE. 0.)       CALL foundBadParameter(32,'Fatal','The fluid conductivity must be positive')
    IF(dx_tentative .LE. 0.)  CALL foundBadParameter(33,'Fatal','dx must be positive')
    IF(TOL_delta_Tg .LE. 0.)  CALL foundBadParameter(34,'Fatal','TOL_delta_Tg must be positive')
    IF(TOL_delta_Tp .LE. 0.)  CALL foundBadParameter(35,'Fatal','TOL_delta_Tp must be positive')
    IF(TOL_delta_Tf .LE. 0.)  CALL foundBadParameter(36,'Fatal','TOL_delta_Tf must be positive')
    IF(TOL_delta_Tb .LE. 0.)  CALL foundBadParameter(37,'Fatal','TOL_delta_Tb must be positive')
      
    IF (ErrorFound()) RETURN
    
    !SET INITIAL VALUES OF STORAGE VARIABLES
    
    Nx = nint(L/dx_tentative)               ! Number of divisions along collector length (x direction), -    
        
    !Setting initial Temperature of collector in storage array
    
    do j = 1,8*(Nx+1)
        call setStaticArrayValue(j,T0)
    end do
    
    i = 1 
    do j = 1,Nx+1 
        call SetStaticArrayValue(i*(Nx+1) + j, T0 + 0.1)
    end do

    i = 5 
    do j = 1,Nx+1 
        call SetStaticArrayValue(i*(Nx+1) + j, T0 + 0.1)
    end do
    
    !  Set the Initial Values of the Outputs
    
    ! this block is to calculate the humidity ratio of the initial inlet (ambient) air based on initial Ta, rh and Patm from weather data
    Ta = getInputValue(3)
    rh = getInputValue(4)
    p_atm = getInputValue(6)
    Tabs = Ta + 273.15
    Pvsat = exp(-5.8002206e3/Tabs + 1.3914493e0 - 4.8640239e-2*Tabs + 4.1764768e-5*Tabs**2 - 1.4452093e-8*Tabs**3 + 6.5459673e0*log(Tabs))
    Pv = (rh/100.)*Pvsat
    Xa = (0.621945*Pv/(p_atm - Pv))
    ! end of this block
    
    Call setOutputValue(1, T0)                  ! the initial outlet air temperature is the initial collector temperature
    Call setOutputValue(2, getInputValue(7))    ! the initial outlet flow rate is the input flow rate
    Call setOutputValue(3, 0.d0)                ! the initial collector efficiency is 0
    Call setOutputValue(4, 0.d0)                ! the useful energy rate
    Call setOutputValue(5, 0.d0)                ! number of iterations
    Call setOutputValue(6, Xa)                  ! initial outlet humidity ratio (from ambient conditions)     
    RETURN 
EndIf

!------------------------------------------------------------------------------ 
! "MULTIPLE UNIT" Manipulations ReRead the Parameters if Another Unit of This Type Has Been Called Last 
If (getIsReReadParameters()) Then 
    !  Read in the Values of the Parameters from the Input File 
    relax       = getParameterValue(1) ! Underrelaxation factor (between 0 and 1)
    L           = getParameterValue(2) ! Collector length (m) 
    W           = getParameterValue(3) ! Collector width (m)
    H           = getParameterValue(4) ! Air channel height (m)
    gap_gp      = getParameterValue(5) ! Gap bewtween cover and absorber (m)
    tilt_deg    = getParameterValue(6) ! Collector tilt from horizontal (deg)
    thick_g     = getParameterValue(7) ! Cover thickness (m)
    thick_p     = getParameterValue(8) ! Absorber thickness (m)
    thick_b     = getParameterValue(9) ! Bottom plate thickness (m)
    k_g         = getParameterValue(10) ! Cover conductivity (W/m-K)
    k_p         = getParameterValue(11) ! Absorber conductivity (W/m-K)
    k_b         = getParameterValue(12) ! Bottom plate conductivity (W/m-K)
    c_g         = getParameterValue(13) ! Cover specific heat (J/kg-K)
    c_p         = getParameterValue(14) ! Absorber specific heat (J/kg-K)
    c_b         = getParameterValue(15) ! Bottom plate specific heat (J/kg-K)
    rho_g       = getParameterValue(16) ! Cover density (kg/m3)
    rho_p       = getParameterValue(17) ! Absorber density (kg/m3)
    rho_b       = getParameterValue(18) ! Bottom plate density (kg/m3)
    thick_i     = getParameterValue(19) ! Insulation thickness (m)
    k_i         = getParameterValue(20) ! Insulation conductivity (W/m-K)
    alpha_g     = getParameterValue(21) ! Cover absorptance (-) 
    alpha_p     = getParameterValue(22) ! Absorber absorptance (-) 
    alpha_b     = getParameterValue(23) ! Bottom plate absorptance (-)
    tau_g       = getParameterValue(24) ! Cover transmittance (-)
    e_g         = getParameterValue(25) ! Cover emittance (-)
    e_p_up      = getParameterValue(26) ! Absorber emittance upwards (-)
    e_p_down    = getParameterValue(27) ! Absorber emittance downwards (-)
    e_b         = getParameterValue(28) ! Bottom plate emittance (-) 
    rho_f       = getParameterValue(29) ! Air density (considered constant) (kg/m3)
    c_f         = getParameterValue(30) ! Air specific heat (considered constant) (J/kg-K)
    mu_f        = getParameterValue(31) ! Air dynamic viscosity (considered constant) (Pa-s)
    k_f         = getParameterValue(32) ! Air thermal conductivity (considered constant) (W/m-K) 
    dx_tentative    = getParameterValue(33) ! tentative (might differ from definitive) grid size along collector length (m) 
    TOL_delta_Tg    = getParameterValue(34) ! tolerance for error between iterations for cover temperature (degC)
    TOL_delta_Tp    = getParameterValue(35) ! tolerance for error between iterations for absorber temperature (degC)
    TOL_delta_Tf    = getParameterValue(36) ! tolerance for error between iterations for fluid temperature (degC)
    TOL_delta_Tb    = getParameterValue(37) ! tolerance for error between iterations for bottom plate temperature (degC)
    T0              = getParameterValue(38) ! Initial collector temperature (degC)
EndIf 

!----------------------------------------------------------------------
! "EVERY TIME STEP" manipulations 

! RETRIEVE STORED VALUES

Nx = nint(L/dx_tentative)               ! Number of divisions along collector length (x direction), -

Tg = 0; Tp = 0; Tf = 0; Tb = 0          ! Setting all elements in vectors to 0. Maybe not needed but just in case
Tg_new = 0; Tp_new = 0; Tf_new = 0; Tb_new = 0          ! Setting all elements in vectors to 0. Maybe not needed but just in case

! Retrieve Tg array
i = 0
do j = 1,Nx+1 
    Tg(j-1) = getStaticArrayValue(i*(Nx+1) + j)
end do   

! Retrieve Tp array
i = 1 
do j = 1,Nx+1 
    Tp(j-1) = getStaticArrayValue(i*(Nx+1) + j)
end do

! Retrieve Tf array
i = 2 
do j = 1,Nx+1 
    Tf(j-1) = getStaticArrayValue(i*(Nx+1) + j)
end do

! Retrieve Tb array
i = 3 
do j = 1,Nx+1 
    Tb(j-1) = getStaticArrayValue(i*(Nx+1) + j)
end do

! Retrieve Tg_new array
i = 4
do j = 1,Nx+1 
    Tg_new(j-1) = getStaticArrayValue(i*(Nx+1) + j)
end do   

! Retrieve Tp_new array
i = 5 
do j = 1,Nx+1 
    Tp_new(j-1) = getStaticArrayValue(i*(Nx+1) + j)
end do

! Retrieve Tf_new array
i = 6 
do j = 1,Nx+1 
    Tf_new(j-1) = getStaticArrayValue(i*(Nx+1) + j)
end do

! Retrieve Tb_new array
i = 7 
do j = 1,Nx+1 
    Tb_new(j-1) = getStaticArrayValue(i*(Nx+1) + j)
end do

!Retrieve Current Inputs to the Model
h_dec =     getInputValue(1)
Irr =	    getInputValue(2)    
Ta =        getInputValue(3)
rh =        getInputValue(4)
v_wind =    getInputValue(5)
p_atm =     getInputValue(6)
mdot_f =    getInputValue(7)
Tf_in =     getInputValue(8)


If ((h_dec < 0) .OR. (h_Dec > 24.0)) Call foundBadInput(1,'Fatal', 'The time in decimal hours must be between 0.0 and 24.0')
If (Irr < 0) Call foundBadInput(2,'Fatal', 'The input solar radiation must be at least 0')
If ((rh < 0) .OR. (rh > 100.))  Call foundBadInput(4,'Fatal', 'The input relative humidity must be between 0 and 100')
If (v_wind < 0) Call foundBadInput(5,'Fatal', 'The input wind speed must be at least 0')
If (p_atm < 0) Call foundBadInput(6,'Fatal', 'The input atmospheric pressure must be at least 0')
If (mdot_f <= 0) Call foundBadInput(7,'Fatal', 'The input flow rate must be positive')  
If (ErrorFound() ) Return

! PERFORM CALCULATIONS

Dh = 4*W*H/(2*(W+H))        ! Hydraulic diameter of flow passage (m)
Ac = L*W                    ! Collector area (m2)
Ub = k_i/thick_i            ! Back loss coefficient (W/m2-K)
th_diff_f = k_f/(c_f*rho_f) ! Air thermal diffusivity (m2/s)  
tilt_rad = tilt_deg*pi/180

do j=0,Nx
    x_position(j) = j*L/Nx   ! Vector with the position of all nodes (m)
end do
    
dx = x_position(1) - x_position(0)             
    
! CALCULATE CONSTANT FACTORS FOR EQUATIONS
Ag = k_g*thick_g/(dx**2); Bg = rho_g*thick_g*c_g
Af = 1./(rho_f*c_f*H)
Ap = k_p*thick_p/(dx**2); Bp = rho_p*thick_p*c_p
Ab = k_b*thick_b/(dx**2); Bb = rho_b*thick_b*c_b
    
hr_pg = 0           ! Initializing vector of radiation heat transfer coefficients between absorber and cover
hr_gs = 0           ! Initializing vector of sky radiation heat transfer coefficients
hr_pb = 0           ! Initializing vector of radiation heat transfer coefficients betweeen absorber and back plate
Ra = 0              ! Initializing vector of Rayleigh number
Nu_nc = 0           ! Initializing vector of natural convection Nusselt number
hnc_pg = 0          ! Initializing vector of natural convection heat transfer coefficient         

delta_Tg = 0        ! Matrix of iteration changes of Tg in each node, initialized to 0, from 0 to Nx_max
delta_Tp = 0.       ! Matrix of iteration changes of Tp in each node, initialized to 0, from 0 to Nx_max
delta_Tf = 0.       ! Matrix of iteration changes of Tf in each node, initialized to 0, from 0 to Nx_max
delta_Tb = 0.       ! Matrix of iteration changes of Tb in each node, initialized to 0, from 0 to Nx_max

do j = 0,Nx         !Set to 1 Values of delta_Tg, delta_Tp, delta_Tf and delta_Tb  where our geometry actually exists (i=0 to Nx)
    delta_Tg(j) = 1.0
    delta_Tp(j) = 1.0
    delta_Tf(j) = 1.0
    delta_Tb(j) = 1.0
end do

iterations = 0
dt = dt_hours*3600

! this block is to calculate the humidity ratio and dew point of the ambient air based on Ta, rh and Patm from weather data
Tabs = Ta + 273.15
Pvsat = exp(-5.8002206e3/Tabs + 1.3914493e0 - 4.8640239e-2*Tabs + 4.1764768e-5*Tabs**2 - 1.4452093e-8*Tabs**3 + 6.5459673e0*log(Tabs))
Pv = (rh/100.)*Pvsat
Xa = (0.621945*Pv/(p_atm - Pv))
Tdp = (6.54 + 14.526*log(Pv/1000) + 0.7389*log(Pv/1000)**2 + 0.09486*log(Pv/1000)**3 + 0.4569*(Pv/1000)**0.1984) !from ASHRAE (2017) Fundamentals pyschrometrics, eq 37
! end of this block
                    
Vdot = mdot_f/rho_f                 ! Airflow rate, m3/s
v = Vdot/(W*H)                      ! Average air speed in air channel (m/s)
spec_Vdot = Vdot/Ac                 ! Airflow rate per unit collector area (m3/s-m2)
Re = v*Dh*rho_f/mu_f                ! Reynolds number for flow in collector
        
if (Re <= 2550.0) then
    Nu = 5.385 + 0.148*Re*H/L
else if (Re <= 10000.0) then
    Nu = 4.4e-4*Re**1.2 + 9.37*Re**0.471*H/L
else
    Nu = 0.03*Re**0.74 + 0.788*Re**0.74*H/L
end if

hc_fb = Nu*k_f/Dh               ! Conv. heat transf. coeff. fluid-back plate (W/m2-K)
hc_pf = Nu*k_f/Dh               ! Conv. heat transf. coeff. absorber-fluid (W/m2-K)
        
Ts = ((Ta + 273.15)*(0.711 + 0.0056*Tdp + 0.000073*Tdp**2 + 0.013*cos(15*h_dec*pi/180))**0.25) - 273.15 ! Sky temperature,  C
hc_gw = 5.7 + 3.8*v_wind        ! Wind conv. heat trans. coeff. (W/m2-K)
        
do while (maxval(delta_Tg) > TOL_delta_Tg .or. maxval(delta_Tp) > TOL_delta_Tp .or. maxval(delta_Tf) > TOL_delta_Tf .or. maxval(delta_Tb) > TOL_delta_Tb)
    
    do j=0,Nx
        hr_pg(j) = sigma*((Tp_new(j)+273)**2 + (Tg_new(j)+273)**2)*((Tp_new(j)+273) + (Tg_new(j)+273))/(1/e_p_up + 1/e_g - 1)
        hr_pb(j) = sigma*((Tp_new(j)+273)**2 + (Tb_new(j)+273)**2)*((Tp_new(j)+273) + (Tb_new(j)+273))/(1/e_p_down + 1/e_b - 1)
        hr_gs(j) = sigma*e_g*((Tg_new(j)+273)**2 + (Ts+273)**2)*((Tg_new(j)+273) + (Ts+273))
        
        Ra(j) = 9.81*(1/((Tp_new(j) + Tg_new(j))/2 + 273))*abs(Tp_new(j) - Tg_new(j))*(gap_gp**3)*rho_f/(mu_f*th_diff_f)
        Nu_nc(j) = 1 + 1.44*(1 - 1708*sin(1.8*tilt_rad)**1.6/(Ra(j)*cos(tilt_rad)))*max(1 - 1708/(Ra(j)*cos(tilt_rad)), 0.0) + max((Ra(j)*cos(tilt_rad)/5830)**(1./3.) - 1, 0.0)  
        hnc_pg(j) = Nu_nc(j)*k_f/gap_gp
    end do
    
    ! LOOP OVER Tg NODES FROM i=0 TO Nx    
    do j = 0,Nx
        if (j>0) then ! This is to be able to use the same equations for internal and boundary node
            jm1 = j-1 
        else
            jm1 = j+1
        end if
                
        if (j<Nx) then 
            jp1 = j+1 
        else 
            jp1 = j-1
        end if
                
        Tg_iter = (Bg*Tg(j) + dt*(Ag*(Tg_new(jm1) + Tg_new(jp1)) + (hnc_pg(j) + hr_pg(j))*Tp_new(j) + alpha_g*Irr + hc_gw*Ta + hr_gs(j)*Ts))/(Bg + dt*(2*Ag + hnc_pg(j) + hc_gw + hr_pg(j) + hr_gs(j)))
            
        Tg_iter = relax*Tg_iter + (1-relax)*Tg_new(j)             ! Applying underrelaxation 
            
        delta_Tg(j) = abs(Tg_iter - Tg_new(j))
                        
        Tg_new(j) = Tg_iter
    end do
            
    ! LOOP OVER Tp NODES FROM i=0 TO Nx
    do j = 0,Nx
                
        if (j>0) then ! This is to be able to use the same equations for internal and boundary node
            jm1 = j-1 
        else
            jm1 = j+1
        end if
                
        if (j<Nx) then 
            jp1 = j+1 
        else 
            jp1 = j-1
        end if
            
        Tp_iter = (Bp*Tp(j) + dt*((hnc_pg(j)+ hr_pg(j))*Tg_new(j) + Ap*(Tp_new(jm1) + Tp_new(jp1)) + hc_pf*Tf_new(j) + hr_pb(j)*Tb_new(j) + alpha_p*tau_g*Irr))/(Bp + dt*(2*Ap + hnc_pg(j) + hc_pf + hr_pg(j) + hr_pb(j)))
    
        Tp_iter = relax*Tp_iter + (1-relax)*Tp_new(j)             ! Applying underrelaxation 
            
        delta_Tp(j) = abs(Tp_iter - Tp_new(j))
            
        Tp_new(j) = Tp_iter
    end do
            
    ! LOOP OVER AIR NODES FROM i=0 TO Nx
    Tf_iter = Tf_in
    Tf_new(0) = Tf_iter
    delta_Tf(0) = Tf_iter - Tf_new(0)
    do j = 1,Nx
        Tf_iter = (Tf(j) + dt*(Af*hc_pf*Tp_new(j) + (v/dx)*Tf_new(j-1) + Af*hc_fb*Tb_new(j)))/(1 + dt*(v/dx + Af*(hc_pf + hc_fb)))
            
        Tf_iter = relax*Tf_iter + (1-relax)*Tf_new(j)             ! Applying underrelaxation 
            
        delta_Tf(j) = abs(Tf_iter - Tf_new(j))
            
        Tf_new(j) = Tf_iter
    end do
            
    ! LOOP OVER Tb NODES FROM i=0 TO Nx
    do j = 0,Nx
            
        if (j>0) then ! This is to be able to use the same equations for internal and boundary node
            jm1 = j-1 
        else
            jm1 = j+1
        end if
                
        if (j<Nx) then 
            jp1 = j+1 
        else 
            jp1 = j-1
        end if 

        Tb_iter = (Bb*Tb(j) + dt*(Ub*Ta + hr_pb(j)*Tp_new(j) + hc_fb*Tf_new(j) + Ab*(Tb_new(jm1) +Tb_new(jp1))))/(Bb + dt*(2*Ab + hr_pb(j) + hc_fb + Ub))
        
        Tb_iter = relax*Tb_iter + (1-relax)*Tb_new(j)             ! Applying underrelaxation 
            
        delta_Tb(j) = abs(Tb_iter - Tb_new(j))
            
        Tb_new(j) = Tb_iter
    end do
        
    iterations = iterations + 1
    
end do            
                            
if (Irr > 0) then
	eff = ((Tf(Nx) - Tf_in)*mdot_f*c_f)/(Irr*Ac)
else
    eff = 0
end if

Usef_Q_rate = mdot_f*c_f*(Tf(Nx) - Tf_in)

! SET STORAGE VALUES

! Set Tg_new array
i = 4
do j = 1,Nx+1 
    call SetStaticArrayValue(i*(Nx+1) + j, Tg_new(j-1))
end do   

! Set Tp_new array
i = 5 
do j = 1,Nx+1 
    call SetStaticArrayValue(i*(Nx+1) + j, Tp_new(j-1))
end do

! Set Tf_new array
i = 6 
do j = 1,Nx+1 
    call SetStaticArrayValue(i*(Nx+1) + j, Tf_new(j-1))
end do

! Set Tb_new array
i = 7 
do j = 1,Nx+1 
    call SetStaticArrayValue(i*(Nx+1) + j, Tb_new(j-1))
end do


! SET OUTPUT VALUES    
Call setOutputValue(1, Tf(Nx))          ! outlet air temperature
Call setOutputValue(2, mdot_f)          ! the flow rate is the input flow rate
Call setOutputValue(3, eff)             ! collector efficiency
Call setOutputValue(4, Usef_Q_rate)     ! useful enery rate (W)
Call setOutputValue(5, iterations)      ! iterations needed in this timestep
Call setOutputValue(6, Xa)          ! outlet air humidity ratio (unchanged from ambient conditions)

Return
End Subroutine TYPE202