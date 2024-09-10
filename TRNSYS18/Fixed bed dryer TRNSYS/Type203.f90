  SUBROUTINE TYPE203

! USE STATEMENTS 
Use TrnsysFunctions 
!-----------------------------------------------------------------------

!export this subroutine for its use in external DLLs. 
!DEC$ATTRIBUTES DLLEXPORT :: TYPE203
    
implicit none

INTERFACE 
   FUNCTION equilibrium_rh(Ta,Xp)
        double precision :: Ta, Xp
   END FUNCTION equilibrium_rh
   
   FUNCTION equilibrium_Xp(Ta,rh)
        double precision :: Ta, rh
   END FUNCTION equilibrium_Xp
   
   FUNCTION drying_rate_K_parameter(Ta,rh,v_a)
        double precision :: Ta, rh, v_a
   END FUNCTION drying_rate_K_parameter
   
   FUNCTION drying_rate_N_parameter()
   END FUNCTION drying_rate_N_parameter   
   
END INTERFACE

! LOCAL VARIABLE DECLARATIONS
    
!   From Trnsys
double precision :: Time, dt_hours ! the current simulation time and the timestep (h)
 

!   What are considered TRNSYS parameters, meaning user inputs that dont change        
double precision :: relax                                                               !1
double precision :: Mass, Xp0, T0, porosity, ap, Dp, c_dp, rho_bulk, Side_length        !9
double precision :: rho_a, c_a, c_v, c_w, mu_a, dz_tentative                            !6
double precision :: TOL_delta_Ta, TOL_delta_Xa, TOL_delta_Tp, TOL_delta_Xp              !4

integer, parameter :: M_max = 500                            ! 1, Maximum number of x nodes along collector length
integer, parameter :: numb_stor_div = 8                ! this is the number of divisions along dryer height to output the variables

double precision :: pi, hfg, Xp0_wetting

!   What are considered TRNSYS inputs, meaning time dependent user inputs
double precision :: h_dec, Ta_in, Xa_in, p_atm, mdot_a  

!   Internal variables necessary which only need to be calculated once
double precision :: rho_dp_bulk, Bulk_volume, Mp0
double precision :: rh0, pvsat0, pv0, Xa0
double precision :: Height, Area_bin, dz, dt

integer :: Static_storage_size, M
double precision :: z_position(0:M_max)

!   Internal variables necessary which may need to be calculated every iteration or timestep
double precision :: Variable_update
double precision :: Ta(0:M_max), Ta_new(0:M_max), Xa(0:M_max), Xa_new(0:M_max), Tp(0:M_max), Tp_new(0:M_max), Xp(0:M_max), Xp_new(0:M_max)
double precision :: delta_Xp(0:M_max), delta_Tp(0:M_max), delta_Xa(0:M_max), delta_Ta(0:M_max)
double precision :: K, N
double precision :: Ta_iter, Xa_iter, Tp_iter, Xp_iter, Xp_eq, Xp_current, Xp_init, MR, teq
double precision :: pvsat_Ta, pv_Ta, pv_Ta_in, pvsat_Ta_in, rh, A, B, C, saturation_vapor_pressure   
double precision :: Vdot_a, v_a, Re, h, Tair_in, mdot_air, Xair_in, rh_in, rho_da, Tabs, iterations
integer :: i, j
integer :: write_pos_index(0:numb_stor_div)

integer :: CurrentUnit, CurrentType                         ! used in Message calls if there is an error (from TRNSYS forum)    
    
! Variables that will be outputs
double precision :: Ta_out     ! mdot_a is also an output but has already been declared as input

Data pi/3.1415927/, hfg/2253430/, Xp0_wetting/0.00/

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
    Mass        = getParameterValue(2) ! Initial product mass at initial moisture content  (kg)
    rho_bulk    = getParameterValue(9) ! Product bulk density at initial moisture (kg/m3)
    Side_length = getParameterValue(10) ! Dryer side length (assuming floor area is square) (m)
    dz_tentative = getParameterValue(16)
    
    Area_bin = Side_length**2    ! Bin floor area, m2 (if square area bin)

    Bulk_volume = Mass/rho_bulk                 ! Volume of bulk, m3 
    
    
    Height = Bulk_volume/Area_bin       ! Height of the grain bed, m
    
   
    M = int(Height/dz_tentative)               ! Number of divisions along bulk height, -    
    
    
    do while (MOD(M,numb_stor_div) .NE. 0)
        M = M + 1
    end do
    
    Static_storage_size = 8*(M+1)           ! Total static storage spots needed for T's and T_new's
    
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
    Mass        = getParameterValue(2) ! Initial product mass at initial moisture content  (kg)
    rho_bulk    = getParameterValue(9) ! Product bulk density at initial moisture (kg/m3)
    Side_length = getParameterValue(10) ! Dryer side length (assuming floor area is square) (m)
    dz_tentative = getParameterValue(16) 
    
    
    Area_bin = Side_length**2    ! Bin floor area, m2 (if square area bin)
    Bulk_volume = Mass/rho_bulk                 ! Volume of bulk, m3 
    
    Height = Bulk_volume/Area_bin       ! Height of the grain bed, m
    
    M = int(Height/dz_tentative)            ! Number of divisions along bulk height, - 
    
    do while (MOD(M,numb_stor_div) .NE. 0)
        M = M + 1
    end do
 
    Static_storage_size = 8*(M+1)           ! Total static storage spots needed for Tp, Xp, Ta, Xa and corresponding new values 


!   Tell the TRNSYS Engine How This Type Works 
    Call SetNumberofParameters(20)            
    Call SetNumberofInputs(4)            
    Call SetNumberofDerivatives(0)            
    Call SetNumberofOutputs(28)            
    Call SetIterationMode(1)     ! 1 indicates the type should be caller every iteration       
    Call SetNumberStoredVariables(Static_storage_size,0) 
    
    !  Set the Correct Input and Output Variable Types 
    
    Call SetInputUnits(1,'TE1')     ! Inlet air temperature (C)
    Call SetInputUnits(2,'DM1')     ! Inlet air humidity ratio (kg/kg)
    Call SetInputUnits(3,'PR3')     ! Atmospheric pressure (Pa)
    Call SetInputUnits(4,'MF2')     ! Mass flow rate (kg/s)
    
    Call SetOutputUnits(1,'DM1')    ! Inlet air relative humidity (-)
    Call SetOutputUnits(2,'TE1')    ! Dryer outlet air temperature (C)    
    Call SetOutputUnits(3,'DM1')    ! Dryer outlet air humidity ratio (kg/kg)
    Call SetOutputUnits(4,'MF2')    ! Mass flow rate (kg/s)
    Call SetOutputUnits(5,'DM1')    ! Avg. product moisture content (kgw/kgdm)
    Call SetOutputUnits(6,'TE1')    ! Avg. product temperature (C)
    Call SetOutputUnits(7,'DM1')    ! Number of iterations
    
    Call SetOutputUnits(8,'DM1')    ! Moisture at 0 height
    Call SetOutputUnits(9,'DM1')    ! Moisture at 0.125 height
    Call SetOutputUnits(10,'DM1')    ! Moisture at 0.25 height
    Call SetOutputUnits(11,'DM1')    ! Moisture at 0.375 height
    Call SetOutputUnits(12,'DM1')    ! Moisture at 0.5 height
    Call SetOutputUnits(13,'DM1')    ! Moisture at 0.625 height
    Call SetOutputUnits(14,'DM1')    ! Moisture at 0.75 height
    Call SetOutputUnits(15,'DM1')    ! Moisture at 0.875 height
    Call SetOutputUnits(16,'DM1')    ! Moisture at 1 height

    Call SetOutputUnits(17,'TE1')    ! Product temperature at 0 height
    Call SetOutputUnits(18,'TE1')    ! Product temperature at 0.125 height
    Call SetOutputUnits(19,'TE1')    ! Product temperature at 0.25 height
    Call SetOutputUnits(20,'TE1')    ! Product temperature at 0.375 height
    Call SetOutputUnits(21,'TE1')    ! Product temperature at 0.5 height
    Call SetOutputUnits(22,'TE1')    ! Product temperature at 0.625 height
    Call SetOutputUnits(23,'TE1')    ! Product temperature at 0.75 height
    Call SetOutputUnits(24,'TE1')    ! Product temperature at 0.875 height
    Call SetOutputUnits(25,'TE1')    ! Product temperature at 1 height
    Call SetOutputUnits(26,'DM1')    ! Bulk height
    Call SetOutputUnits(27,'DM1')    ! M
    Call SetOutputUnits(28,'DM1')    ! dz
      
    Return 
EndIf 
    
    
 !------------------------------------------------------------------------------ 
    !DO ALL OF THE "START TIME" MANIPULATIONS HERE - There Are No Iterations at the Intial Time 
If (getIsStartTime()) Then
    !  Read in the Values of the Parameters from the Input File 
    relax       = getParameterValue(1) ! Underrelaxation factor (between 0 and 1)
    Mass        = getParameterValue(2) ! Initial product mass at initial moisture content  (kg)
    Xp0         = getParameterValue(3) ! Initial moisture content in dry basis (-)
    T0          = getParameterValue(4) ! Initial temperature (C)
    porosity    = getParameterValue(5) ! Bulk porosity (-)
    ap          = getParameterValue(6) ! Surface area of product per unit bulk volume (m2/m3) 
    Dp          = getParameterValue(7) ! Equivalent particle diameter (m)
    c_dp        = getParameterValue(8) ! Product specific heat in dry basis (J/kg-K)
    rho_bulk    = getParameterValue(9) ! Product bulk density at initial moisture (kg/m3)
    Side_length = getParameterValue(10) ! Dryer side length (assuming floor area is square) (m)
    rho_a       = getParameterValue(11) ! Air density (considered constant) (kg/m3)
    c_a         = getParameterValue(12) ! Air specific heat (considered constant) (J/kg-K)
    mu_a        = getParameterValue(13) ! Air dynamic viscosity (considered constant) (Pa-s)
    c_v         = getParameterValue(14) ! Specific heat of water vapor (J/kg-K) 
    c_w         = getParameterValue(15) ! Specific heat of liquid water (J/kg-K) 
    dz_tentative    = getParameterValue(16) ! tentative (might differ from definitive) grid size along dryer height (m)
    TOL_delta_Xp    = getParameterValue(17) ! tolerance for error between iterations for product moisture (kg/kg)
    TOL_delta_Tp    = getParameterValue(18) ! tolerance for error between iterations for product (C)
    TOL_delta_Xa    = getParameterValue(19) ! tolerance for error between iterations for air humidity ratio (kg/kg)
    TOL_delta_Ta    = getParameterValue(20) ! tolerance for error between iterations for air temperature (C)
    
   !   Check the parameters for problems and RETURN if any are found 
    IF((relax .GT. 1.).OR.(relax .LE. 0.)) CALL foundBadParameter(1,'Fatal','The underrelax factor must be above 0 and up to 1')
    IF(Mass .LE. 0.)      CALL foundBadParameter(2,'Fatal','The initial product mass must more than 0') 
    IF(Xp0 .LE. 0.)       CALL foundBadParameter(3,'Fatal','The initial product moisture content must be more than 0')
    IF((porosity .LE. 0.) .OR. (porosity .LE. 0))   CALL foundBadParameter(5,'Fatal','The product bulk porosity must be above 0 and up to 1')
    IF(ap .LE. 0.)        CALL foundBadParameter(6,'Fatal','The bulk specific surface area must be more than 0')
    IF(Dp .LE. 0.)        CALL foundBadParameter(7,'Fatal','The product equivalent particle diameter must be more than 0')
    IF(c_dp .LE. 0.)       CALL foundBadParameter(8,'Fatal','The product specific heat must be more than 0')
    IF(rho_bulk .LE. 0.)       CALL foundBadParameter(9,'Fatal','The product bulk density must be more than 0')        
    IF(Side_length .LE. 0.)       CALL foundBadParameter(10,'Fatal','The side length of the dryer must be more than 0')
    IF(rho_a .LE. 0.)       CALL foundBadParameter(11,'Fatal','The air density must be more than 0')
    IF(c_a .LE. 0.)       CALL foundBadParameter(12,'Fatal','The air specific heat must be more than 0')
    IF(mu_a .LE. 0.)     CALL foundBadParameter(13,'Fatal','The air viscosity must be more than 0')
    IF(c_v .LE. 0.)     CALL foundBadParameter(14,'Fatal','The water vapor specific heat must be more than 0')
    IF(c_w .LE. 0.)     CALL foundBadParameter(15,'Fatal','The liquid water specific heat must be more than 0')
    IF(dz_tentative .LE. 0.)    CALL foundBadParameter(16,'Fatal','dz_tentative must be more than 0 ')
    IF(TOL_delta_Xp .LE. 0.)    CALL foundBadParameter(17,'Fatal','TOL_delta_Xp must be more than 0')
    IF(TOL_delta_Tp .LE. 0.)    CALL foundBadParameter(18,'Fatal','TOL_delta_Tp must be more than 0')
    IF(TOL_delta_Xa .LE. 0.)    CALL foundBadParameter(19,'Fatal','TOL_delta_Xa must be more than 0')
    IF(TOL_delta_Ta .LE. 0.)    CALL foundBadParameter(20,'Fatal','TOL_delta_Ta must be more than 0')
      
    IF (ErrorFound()) RETURN
    
    !SET INITIAL VALUES OF STORAGE VARIABLES
    Area_bin = Side_length**2    ! Bin floor area, m2 (if square area bin)
    Bulk_volume = Mass/rho_bulk                 ! Volume of bulk, m3 
    Height = Bulk_volume/Area_bin       ! Height of the grain bed, m
    M = nint(Height/dz_tentative)            ! Number of divisions along bulk height, - 
    
    p_atm = getInputValue(3)
    
    do while (MOD(M,numb_stor_div) .NE. 0)
        M = M + 1
    end do
    
    rh0 = equilibrium_rh(T0,Xp0)           ! Relative humidity of intersticial air in equilibrium with product, -
    
    
    Tabs = T0 + 273.15
    
    pvsat0 = exp(-5.8002206e3/Tabs + 1.3914493e0 - 4.8640239e-2*Tabs + 4.1764768e-5*Tabs**2 - 1.4452093e-8*Tabs**3 + 6.5459673e0*log(Tabs))
    pv0 = rh0*pvsat0                            ! Actual initial vapor pressure of the intersticial air, Pa
    
    
    Xa0 = 0.621945*pv0/(p_atm - pv0)                ! Initial humidity ratio of intersticial air, kg/kg
        
    !Setting initial temperature (product and air), moisture and humidity in dryer, in that order
    i=0
    do j = 1,M+1
        call setStaticArrayValue(i*(M+1) + j,Xp0)
    end do
    
    i=1
    do j = 1,M+1
        call setStaticArrayValue(i*(M+1) + j,T0)
    end do
    
    i=2
    do j = 1,M+1
        call setStaticArrayValue(i*(M+1) + j,Xa0)
    end do
    
    i=3
    do j = 1,M+1
        call setStaticArrayValue(i*(M+1) + j,T0)
    end do
    
    i=4
    do j = 1,M+1
        call setStaticArrayValue(i*(M+1) + j,Xp0)
    end do
    
    i=5
    do j = 1,M+1
        call setStaticArrayValue(i*(M+1) + j,T0)
    end do
    
    i=6
    do j = 1,M+1
        call setStaticArrayValue(i*(M+1) + j,Xa0)
    end do
    
    i=7
    do j = 1,M+1
        call setStaticArrayValue(i*(M+1) + j,T0)
    end do
       
    !  Set the Initial Values of the Outputs  
    Call setOutputValue(1, rh0)             ! The inlet rh at time 0 
    Call setOutputValue(2, T0)              ! The initial outlet air temperature is T0
    Call setOutputValue(3, Xa0)             ! The initial outlet humidity is Xa0
    Call setOutputValue(4, getInputValue(4))! The initial mass flow rate is the initial input value
    Call setOutputValue(5, Xp0)             ! The initial average moisture content is Xp0
    Call setOutputValue(6, T0)              ! The initial average temperature is T0
    Call SetOutputValue(7, 0.d0)            ! The initial number of iterations is 0
    
    Call SetOutputValue(8, Xp0)    ! Moisture at 0 height
    Call SetOutputValue(9, Xp0)    ! Moisture at 0.125 height
    Call SetOutputValue(10, Xp0)    ! Moisture at 0.25 height
    Call SetOutputValue(11, Xp0)    ! Moisture at 0.375 height
    Call SetOutputValue(12, Xp0)    ! Moisture at 0.5 height
    Call SetOutputValue(13, Xp0)    ! Moisture at 0.625 height
    Call SetOutputValue(14, Xp0)    ! Moisture at 0.75 height
    Call SetOutputValue(15, Xp0)    ! Moisture at 0.875 height
    Call SetOutputValue(16, Xp0)    ! Moisture at 1 height

    Call SetOutputValue(17, T0)    ! Product temperature at 0 height
    Call SetOutputValue(18, T0)    ! Product temperature at 0.125 height
    Call SetOutputValue(19, T0)    ! Product temperature at 0.25 height
    Call SetOutputValue(20, T0)    ! Product temperature at 0.375 height
    Call SetOutputValue(21, T0)    ! Product temperature at 0.5 height
    Call SetOutputValue(22, T0)    ! Product temperature at 0.625 height
    Call SetOutputValue(23, T0)    ! Product temperature at 0.75 height
    Call SetOutputValue(24, T0)    ! Product temperature at 0.875 height
    Call SetOutputValue(25, T0)    ! Product temperature at 1 height
    
    Call SetOutputValue(26, 0)    ! Height
    Call SetOutputValue(27, 10)    ! M
    Call SetOutputValue(28, 0.008)    ! dz
    
    RETURN 
EndIf

!------------------------------------------------------------------------------ 
! "MULTIPLE UNIT" Manipulations ReRead the Parameters if Another Unit of This Type Has Been Called Last 
If (getIsReReadParameters()) Then 
    !  Read in the Values of the Parameters from the Input File 
    relax       = getParameterValue(1) ! Underrelaxation factor (between 0 and 1)
    Mass        = getParameterValue(2) ! Initial product mass at initial moisture content  (kg)
    Xp0         = getParameterValue(3) ! Initial moisture content in dry basis (-)
    T0          = getParameterValue(4) ! Initial temperature (°C)
    porosity    = getParameterValue(5) ! Bulk porosity (-)
    ap          = getParameterValue(6) ! Surface area of product per unit bulk volume (m2/m3) 
    Dp          = getParameterValue(7) ! Equivalent particle diameter (m)
    c_dp        = getParameterValue(8) ! Product specific heat in dry basis (J/kg-K)
    rho_bulk    = getParameterValue(9) ! Product bulk density at initial moisture (kg/m3)
    Side_length = getParameterValue(10) ! Dryer side length (assuming floor area is square) (m)
    rho_a       = getParameterValue(11) ! Air density (considered constant) (kg/m3)
    c_a         = getParameterValue(12) ! Air specific heat (considered constant) (J/kg-K)
    mu_a        = getParameterValue(13) ! Air dynamic viscosity (considered constant) (Pa-s)
    c_v         = getParameterValue(14) ! Specific heat of water vapor (J/kg-K) 
    c_w         = getParameterValue(15) ! Specific heat of liquid water (J/kg-K) 
    dz_tentative    = getParameterValue(16) ! tentative (might differ from definitive) grid size along dryer height (m)
    TOL_delta_Xp    = getParameterValue(17) ! tolerance for error between iterations for product moisture (kg/kg)
    TOL_delta_Tp    = getParameterValue(18) ! tolerance for error between iterations for product (C)
    TOL_delta_Xa    = getParameterValue(19) ! tolerance for error between iterations for air humidity ratio (kg/kg)
    TOL_delta_Ta    = getParameterValue(20) ! tolerance for error between iterations for air temperature (C)
    
EndIf 

!----------------------------------------------------------------------
! "EVERY TIME STEP" manipulations 

! RETRIEVE STORED VALUES

Area_bin = Side_length**2    ! Bin floor area, m2 (if square area bin)
Bulk_volume = Mass/rho_bulk                 ! Volume of bulk, m3 
Height = Bulk_volume/Area_bin       ! Height of the grain bed, m
M = nint(Height/dz_tentative)               ! Number of divisions along collector length (x direction), -

do while (MOD(M,numb_stor_div) .NE. 0)
    M = M + 1
end do

Ta = 0; Tp = 0; Xa = 0; Xp = 0          ! Setting all elements in vectors to 0. Maybe not needed but just in case
Ta_new = 0; Tp_new = 0; Xa_new = 0; Xp_new = 0          ! Setting all elements in vectors to 0. Maybe not needed but just in case

! Retrieve Xp array
i = 0
do j = 1,M+1 
    Xp(j-1) = getStaticArrayValue(i*(M+1) + j)
end do   

! Retrieve Tp array
i = 1 
do j = 1,M+1 
    Tp(j-1) = getStaticArrayValue(i*(M+1) + j)
end do

! Retrieve Xa array
i = 2 
do j = 1,M+1 
    Xa(j-1) = getStaticArrayValue(i*(M+1) + j)
end do

! Retrieve Ta array
i = 3 
do j = 1,M+1 
    Ta(j-1) = getStaticArrayValue(i*(M+1) + j)
end do

! Retrieve Xp_new array
i = 4
do j = 1,M+1 
    Xp_new(j-1) = getStaticArrayValue(i*(M+1) + j)
end do   

! Retrieve Tp_new array
i = 5 
do j = 1,M+1 
    Tp_new(j-1) = getStaticArrayValue(i*(M+1) + j)
end do

! Retrieve Xa_new array
i = 6 
do j = 1,M+1 
    Xa_new(j-1) = getStaticArrayValue(i*(M+1) + j)
end do

! Retrieve Ta_new array
i = 7 
do j = 1,M+1 
    Ta_new(j-1) = getStaticArrayValue(i*(M+1) + j)
end do

!Retrieve Current Inputs to the Model
Ta_in =     getInputValue(1)
Xa_in =	    getInputValue(2)    
p_atm =     getInputValue(3)
mdot_air =  getInputValue(4)

If (Xa_in <= 0)  Call foundBadInput(2,'Fatal', 'The inlet humidity ratio must positive')
If (p_atm < 0) Call foundBadInput(3,'Fatal', 'The input atmospheric pressure must positive')
If (mdot_air <= 0) Call foundBadInput(4,'Fatal', 'The input flow rate must be positive')    
If (ErrorFound() ) Return


Tabs = Ta_in + 273.15


pvsat_Ta = exp(-5.8002206e3/Tabs + 1.3914493e0 - 4.8640239e-2*Tabs + 4.1764768e-5*Tabs**2 - 1.4452093e-8*Tabs**3 + 6.5459673e0*log(Tabs))    ! Saturation vapor 
                                                                                                    ! pressure at the air temperature in the current node (Pa)                              


pv_Ta = Xa_in*p_atm/(0.621945 + Xa_in)          ! Actual vapor pressure of air in current node (Pa)


rh_in = pv_Ta/pvsat_Ta 

! PERFORM CALCULATIONS

Mp0 = Xp0/(1+Xp0)                   ! Initial product moisture content, wb, kg/kg
rho_dp_bulk = rho_bulk*(1-Mp0)      ! Product bulk dry matter density, kg/m3
                                    ! Since the volume of the product bulk is fixed in the model the dry matter bulk density is constant

dz = Height/M              ! Mesh size (actual), m

do i=0,numb_stor_div
    write_pos_index(i) = int(i*1.d0/(numb_stor_div)*M)
end do
    
! CALCULATE CONSTANT FACTORS FOR EQUATIONS

delta_Xp = 0.        ! Matrix of iteration changes of Xp in each node, initialized to 0, from 0 to M_max
delta_Tp = 0.       ! Matrix of iteration changes of Tp in each node, initialized to 0, from 0 to M_max
delta_Ta = 0.       ! Matrix of iteration changes of Ta in each node, initialized to 0, from 0 to M_max
delta_Xa = 0.       ! Matrix of iteration changes of Xa in each node, initialized to 0, from 0 to M_max

do j = 0,M         !Set to 1 Values of delta_Xp, delta_Tp, delta_Ta and delta_Xa where there is grid (from 0 to M)
    delta_Xp(j) = 1.0
    delta_Tp(j) = 1.0
    delta_Ta(j) = 1.0
    delta_Xa(j) = 1.0
end do

iterations = 0
dt = dt_hours*3600



Vdot_a = mdot_air/rho_a     ! Airflow rate, m3/s
v_a = Vdot_a/Area_bin       ! Superficial velocity of drying air, m/s


Re = v_a*rho_a*Dp/mu_a      ! Reynolds number

h = 0.2755*(c_a/1000.)*(v_a*3600)*rho_a*Re**(-0.34)  ! Heat transfer coefficient according to Brooker, B-A book and paper de B-A 1967


A = dt*v_a/(porosity*dz)    ! Repeated factor in eq for Ta, Xa
             
do while (maxval(delta_Xp) > TOL_delta_Xp .or. maxval(delta_Tp) > TOL_delta_Tp .or. maxval(delta_Xa) > TOL_delta_Xa .or. maxval(delta_Ta) > TOL_delta_Ta)
    ! LOOP OVER Ta NODES FROM j=0 TO M
    Ta_iter = Ta_in
    Ta_new(0) = Ta_iter
    delta_Ta(0) = Ta_iter - Ta_new(0)
    
    do j = 1,M
        if  ((0.287042*(Ta_new(j) + 273.15)*(1 + 1.607858*Xa_new(j))/(0.001*p_atm)) == 0) CALL foundBadParameter(23,'Fatal',' (0.287042*(Ta_new(j) + 273.15)*(1 + 1.607858*Xa_new(j))/(0.001*p_atm)) equal to zero')
        rho_da = 1./(0.287042*(Ta_new(j) + 273.15)*(1 + 1.607858*Xa_new(j))/(0.001*p_atm)) ! Calculate dry air density for current node. First trying with 
                                                                                            !Ta[j] and Xa[j] instead of new values.
                            
        B = porosity*rho_da*(c_a + c_v*Xa_new(j)) ! Repeated denominator of eq for Ta
        
        Ta_iter = (Ta(j) + A*Ta_new(j-1) + dt*h*ap*Tp_new(j)/B)/(1. + A + dt*h*ap/B)   
                                                                                                                                                                      
        Ta_iter = relax*Ta_iter + (1-relax)*Ta_new(j)             ! Applying underrelaxation 
        delta_Ta(j) = abs(Ta_iter - Ta_new(j))
        Ta_new(j) = Ta_iter            
    end do
       
    ! LOOP OVER Xa NODES FROM j=0 TO M
    Xa_iter = Xa_in
    Xa_new(0) = Xa_iter
    delta_Xa(0) = Xa_iter - Xa_new(0)
        
    do j = 1,M
        if  ((0.287042*(Ta_new(j) + 273.15)*(1 + 1.607858*Xa_new(j))/(0.001*p_atm)) == 0) CALL foundBadParameter(23,'Fatal',' (0.287042*(Ta_new(j) + 273.15)*(1 + 1.607858*Xa_new(j))/(0.001*p_atm)) equal to zero')
        rho_da = 1/(0.287042*(Ta_new(j) + 273.15)*(1 + 1.607858*Xa_new(j))/(0.001*p_atm)) ! Calculate dry air density for current node
        
        Xa_iter = (Xa(j) + A*Xa_new(j-1) - rho_dp_bulk/(porosity*rho_da)*(Xp_new(j) - Xp(j)))/(1 + A)   ! Here also using Xp upwind values since air humidity will be
                                                                                            ! more influenced by the product moisture content upwind of the air
                                                                                               ! node than by the downwind product node
        Xa_iter = relax*Xa_iter + (1-relax)*Xa_new(j)             ! Applying underrelaxation 
        delta_Xa(j) = abs(Xa_iter - Xa_new(j))
        Xa_new(j) = Xa_iter
    end do
            
    ! LOOP OVER Tp NODES FROM j=0 TO M
    j = 0
    
    rho_da = 1/(0.287042*(Ta_new(j) + 273.15)*(1 + 1.607858*Xa_new(j))/(0.001*p_atm)) ! Calculate dry air density for current node
    
    
   
    C = (1 + (h*ap*dt*dz - rho_da*v_a*c_v*dt*(Xa_new(j+1) - Xa_new(j)))/(rho_dp_bulk*dz*(c_dp + c_w*Xp_new(j))))
          
    
    Tp_iter = (Tp(j) + (dt/(rho_dp_bulk*(c_dp + c_w*Xp_new(j))))*(h*ap*Ta_new(j) - (hfg + c_v*Ta_new(j))*(rho_da*v_a*(Xa_new(j+1) - Xa_new(j))/dz)))/C
               
    Tp_iter = relax*Tp_iter + (1-relax)*Tp_new(j)             ! Applying underrelaxation  
    delta_Tp(j) = abs(Tp_iter - Tp_new(j))
    Tp_new(j) = Tp_iter
            
    do j = 1,M-1
            
        rho_da = 1/(0.287042*(Ta_new(j) + 273.15)*(1 + 1.607858*Xa_new(j))/(0.001*p_atm)) ! Calculate dry air density for current node
        
        
       
        C = (1 + (2*h*ap*dt*dz - rho_da*v_a*c_v*dt*(Xa_new(j+1) - Xa_new(j-1)))/(2*rho_dp_bulk*dz*(c_dp + c_w*Xp_new(j))))
             
        
      
        Tp_iter = (Tp(j) + (dt/(rho_dp_bulk*(c_dp + c_w*Xp_new(j))))*(h*ap*Ta_new(j) - (hfg + c_v*Ta_new(j))*(rho_da*v_a*(Xa_new(j+1) - Xa_new(j-1))/(2*dz))))/C
                
        Tp_iter = relax*Tp_iter + (1-relax)*Tp_new(j)             ! Applying underrelaxation    
        delta_Tp(j) = abs(Tp_iter - Tp_new(j))
        Tp_new(j) = Tp_iter
    end do
            
    j = M
    
    rho_da = 1/(0.287042*(Ta_new(j) + 273.15)*(1 + 1.607858*Xa_new(j))/(0.001*p_atm)) ! Calculate dry air density for current node
    
   
    C = (1 + (h*ap*dt*dz - rho_da*v_a*c_v*dt*(Xa_new(j) - Xa_new(j-1)))/(rho_dp_bulk*dz*(c_dp + c_w*Xp_new(j))))
    
   
    Tp_iter = (Tp(j) + (dt/(rho_dp_bulk*(c_dp + c_w*Xp_new(j))))*(h*ap*Ta_new(j) - (hfg + c_v*Ta_new(j))*(rho_da*v_a*(Xa_new(j) - Xa_new(j-1))/dz)))/C
        
    Tp_iter = relax*Tp_iter + (1-relax)*Tp_new(j)             ! Applying underrelaxation
    delta_Tp(j) = abs(Tp_iter - Tp_new(j))
    Tp_new(j) = Tp_iter
    
    ! LOOP OVER PRODUCT MOISTURE CONTENT FROM i=0 TO M        
     do j = 0,M
        Tabs = Ta_new(j) + 273.15
        
        pvsat_Ta = exp(-5.8002206e3/Tabs + 1.3914493e0 - 4.8640239e-2*Tabs + 4.1764768e-5*Tabs**2 - 1.4452093e-8*Tabs**3 + 6.5459673e0*log(Tabs))    ! Saturation 
                                                                                                    !vapor pressure at the air temperature in the current node (Pa)
                                       
        pv_Ta = Xa_new(j)*p_atm/(0.621945 + Xa_new(j))          ! Actual vapor pressure of air in current node (Pa)
                
        rh = pv_Ta/pvsat_Ta 
                              
        if (rh >= 0.99) then
            rh = 0.99
        end if
               
        Xp_eq = equilibrium_Xp(Ta_new(j),rh) ! Equilibruim moisture content
                       
        K = drying_rate_K_parameter(Ta_new(j),rh,v_a)   ! Thin layer Page drying K constant for wheat from Ramaj (2021)
        N = drying_rate_N_parameter()                   ! Thin layer Page drying N constant for wheat from Ramaj (2021)
            
        Xp_current = Xp_new(j)
        Xp_init = Xp0
            
        if (Xp_current > Xp_eq) then          ! In this case drying should happen in that layer. We assume that Xp_current is never exactly the X_eq
            if (Xp_eq == Xp_init) then
                Xp_eq = Xp_init - 0.0001
            end if
                    
            if (Xp_current > Xp_init) then
                Xp_init = Xp_current
            end if

            MR = (Xp_current-Xp_eq)/(Xp_init - Xp_eq)
                    
            teq = (-log(MR)/K)**(1./N)*60.
                       
        else
            MR = (Xp_current-Xp_eq)/(Xp0_wetting - Xp_eq)    ! In case Xp_current is less than Xp_eq, wetting should happen, and we use a special initial moisture content
                  
            teq = (-log(MR)/K)**(1./N)*60.
        end if
        
        Xp_iter = (Xp(j) + (dt/60.)*N*K*((teq+dt)/60.)**(N-1)*Xp_eq)/(1 + (dt/60.)*N*K*((teq+dt)/60.)**(N-1))
                
        Xp_iter = relax*Xp_iter + (1-relax)*Xp_new(j)             ! Applying underrelaxation  
            
        delta_Xp(j) = abs(Xp_iter - Xp_new(j))
            
        Xp_new(j) = Xp_iter
            
     end do
    
     iterations = iterations + 1

end do


! SET STORAGE VALUES

! Set Xp_new array
i = 4
do j = 1,M+1 
    call SetStaticArrayValue(i*(M+1) + j, Xp_new(j-1))
end do   

! Set Tp_new array
i = 5 
do j = 1,M+1 
    call SetStaticArrayValue(i*(M+1) + j, Tp_new(j-1))
end do

! Set Xa_new array
i = 6 
do j = 1,M+1 
    call SetStaticArrayValue(i*(M+1) + j, Xa_new(j-1))
end do

! Set Ta_new array
i = 7 
do j = 1,M+1 
    call SetStaticArrayValue(i*(M+1) + j, Ta_new(j-1))
end do


! SET OUTPUT VALUES    
Call setOutputValue(1, rh_in)           ! The inlet rh at current time step 
Call setOutputValue(2, Ta_new(M))       ! The outlet air temperature
Call setOutputValue(3, Xa_new(M))        ! The outlet humidity ratio
Call setOutputValue(4, mdot_air)        ! The mass flow rate at current time step
Call setOutputValue(5, sum(Xp_new)/(M+1)) ! The average moisture content
Call setOutputValue(6, sum(Tp_new)/(M+1)) ! The average product temperature
Call SetOutputValue(7, iterations)      ! The initial number of iterations is 0
    
Call SetOutputValue(8, Xp_new(write_pos_index(0)))    ! Moisture at 0 height
Call SetOutputValue(9, Xp_new(write_pos_index(1)))    ! Moisture at 0.125 height
Call SetOutputValue(10, Xp_new(write_pos_index(2)))    ! Moisture at 0.25 height
Call SetOutputValue(11, Xp_new(write_pos_index(3)))    ! Moisture at 0.375 height
Call SetOutputValue(12, Xp_new(write_pos_index(4)))    ! Moisture at 0.5 height
Call SetOutputValue(13, Xp_new(write_pos_index(5)))    ! Moisture at 0.625 height
Call SetOutputValue(14, Xp_new(write_pos_index(6)))    ! Moisture at 0.75 height
Call SetOutputValue(15, Xp_new(write_pos_index(7)))    ! Moisture at 0.875 height
Call SetOutputValue(16, Xp_new(write_pos_index(8)))    ! Moisture at 1 height

Call SetOutputValue(17, Tp_new(write_pos_index(0)))    ! Product temperature at 0 height
Call SetOutputValue(18, Tp_new(write_pos_index(1)))    ! Product temperature at 0.125 height
Call SetOutputValue(19, Tp_new(write_pos_index(2)))    ! Product temperature at 0.25 height
Call SetOutputValue(20, Tp_new(write_pos_index(3)))    ! Product temperature at 0.375 height
Call SetOutputValue(21, Tp_new(write_pos_index(4)))    ! Product temperature at 0.5 height
Call SetOutputValue(22, Tp_new(write_pos_index(5)))    ! Product temperature at 0.625 height
Call SetOutputValue(23, Tp_new(write_pos_index(6)))    ! Product temperature at 0.75 height
Call SetOutputValue(24, Tp_new(write_pos_index(7)))    ! Product temperature at 0.875 height
Call SetOutputValue(25, Tp_new(write_pos_index(8)))    ! Product temperature at 1 height

Call SetOutputValue(26, Height)    ! Height
Call SetOutputValue(27, real(M))    ! M
Call SetOutputValue(28, dz)    ! dz

    

Return
End Subroutine TYPE203


! FUNCTION TO CALCULATE THE EQUILIBRIUM RELATIVE HUMIDITY BASED ON SORPTION ISOTHERM EQUATION
double precision function equilibrium_rh(Ta,Xp)
    implicit none
    double precision :: Ta, Xp
    equilibrium_rh = (1./(((0.129 - 6.46e-4*Ta)/Xp)**2.944 + 1.))
end function equilibrium_rh

! FUNCTION TO CALCULATE THE EQUILIBRIUM MOISTURE CONTENT BASED ON SORPTION ISOTHERM EQUATION
double precision function equilibrium_Xp(Ta,rh)
    implicit none
    double precision :: Ta, rh
    equilibrium_Xp = (0.129 - 6.46e-4*Ta)*(rh/(1.-rh))**(1./2.944)
end function equilibrium_Xp 
    
! FUNCTION TO CALCULATE THE K DRYING PARAMETER OF THE THIN LAYER DRYING EQUATION (HERE FOR THE PAGE MODEL)    
double precision function drying_rate_K_parameter(Ta,rh,v_a)
    implicit none
    double precision :: Ta, rh, v_a
    drying_rate_K_parameter = 2.8e-3*exp(0.059*Ta)*(100*rh)**(-0.139)*v_a**0.025     ! Page k parameter for wheat from Ramaj, 2021
end function drying_rate_K_parameter
    
! FUNCTION TO CALCULATE THE N DRYING PARAMETER OF THE THIN LAYER DRYING EQUATION (HERE FOR THE PAGE MODEL)
! IN THIS CASE N IS A CONSTANT BUT OFTEN IT IS NOT    
double precision function drying_rate_N_parameter()
    implicit none
    drying_rate_N_parameter = 0.784     ! Page N parameter for wheat from Ramaj, 2021
end function drying_rate_N_parameter