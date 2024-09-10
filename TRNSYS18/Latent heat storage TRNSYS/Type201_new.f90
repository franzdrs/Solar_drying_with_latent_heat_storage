SUBROUTINE TYPE201

! USE STATEMENTS 
Use TrnsysFunctions 
!-----------------------------------------------------------------------

!export this subroutine for its use in external DLLs. 
!DEC$ATTRIBUTES DLLEXPORT :: TYPE201
    
implicit none
    
! LOCAL VARIABLE DECLARATIONS
    
!   From Trnsys
double precision :: Time, dt_hours ! the current simulation time and the timestep (h)
    
!   What are considered TRNSYS parameters, meaning user inputs that dont change
double precision :: relax, rho_s, rho_l, ks, kl, cs, cl, L, Tm, T_range   ! 10 
double precision :: rho_a, c_a, mu_a, k_a, Pr                             ! 5
double precision :: rho_w, thickn_w, c_w, k_w                             ! 4
double precision :: Slab_length, Slab_width_s, Slab_thickness, Air_gap    ! 4
double precision :: T0                                                    ! 1
double precision :: dx_tentative, dy_tentative                            ! 2
integer :: N_slabs                                                        ! 1
double precision :: TOL_delta_H, TOL_delta_Ta, TOL_delta_Tw               ! 3
integer, parameter :: M_max = 100, N_max = 10                             ! 2, Maximum number of x and y nodes for air, wall, slab

double precision :: H0        ! this is not a parameter to be input by user, just the initial enthalpy calculated from T0

!double precision :: Store_dt_min, Store_dt_min_details
    
!   What are considered TRNSYS inputs, meaning time dependent user inputs
double precision :: Tair_in, mdot_a, Xa_in
    
!   Internal variables necessary which only need to be calculated once
double precision :: expansion, Mass_pcm, Hs, Hl, Ts, Tl, A, fe, Gap_AR          ! 9
double precision :: Slab_half_thickness, Half_air_gap, Dh, Duct_width           ! 4
integer :: N_channels, M, N, Static_storage_size                                ! 4
double precision :: dx, dy, dt                                                      ! 3

double precision :: x_position(0:M_max), y_position(0:N_max), grid_weights(0:N_max,0:M_max)  ! 3


!   Internal variables necessary which may need to be calculated every iteration or timestep 
double precision :: Re, Nu, ff, Xaster, h_aw, W, R, F, pcm_width, k_int, va, S, V_a
double precision :: kx_ant, kx_post, ky_ant, ky_post
double precision :: Ta_iter, Tw_iter, T_iter, H_iter, Variable_update, iterations
integer :: i, j, jm1, jp1 !, iterations
double precision :: Ta(0:M_max), Ta_new(0:M_max) 
double precision :: Tw(0:M_max), Tw_new(0:M_max)
double precision :: T(0:N_max,0:M_max), T_new(0:N_max,0:M_max) 
double precision :: H(0:N_max,0:M_max), H_new(0:N_max,0:M_max), Lambda(0:N_max,0:M_max)
 
double precision :: k(0:N_max,0:M_max)
double precision :: delta_H(0:N_max,0:M_max), delta_Ta(0:M_max), delta_Tw(0:M_max)
integer :: CurrentUnit, CurrentType                         ! used in Message calls if there is an error (from TRNSYS forum)
    
! Variables that will be outputs
double precision :: Ta_out_mixed, Lambda_avg, Total_enthalpy, Avg_temperature    ! mdot_a is also an output but has already been declared as input
       
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
    
    Slab_length     = getParameterValue(21) ! Slab length, m
    Slab_thickness  = getParameterValue(23) ! Slab thickness, m
    dx_tentative    = getParameterValue(25) ! tentative (might differ from definitive) grid size along slab length , m
    dy_tentative    = getParameterValue(26) ! tentative (might differ from definitive) grid size along slab thickness, m
    Slab_half_thickness = Slab_thickness/2.0 ! Slab half thickness, m
    
    M = nint(Slab_length/dx_tentative)              ! Number of divisions along slab length (x direction), -   
    N = nint(Slab_half_thickness/dy_tentative)      ! Number of divisions along half the slab thickness (y direction), - 
    
    Static_storage_size = 4*(M+1) + 4*(M+1)*(N+1)   ! Total static storage spots needed for Ta, Ta_new, Tw, Tw_new, H, H_new, T and T_new
    
    !   UPDATE STORED VALUES FROM Ta_new, Tw_new, T_new and H_new to Ta, Tw, T, and H
    
    do j = Static_storage_size/2+1, Static_storage_size 
        Variable_update = getStaticArrayValue(j)
        call setStaticArrayValue(j-Static_storage_size/2, Variable_update)
    end do
        
    Return 
Endif

!---------------------------------------------------------------------- 
!   DO ALL OF THE "INITIALIZATION CALL  MANIPULATIONS HERE 
If (getIsFirstCallofSimulation()) Then
    
    Slab_length     = getParameterValue(21) ! Slab length, m
    Slab_thickness  = getParameterValue(23) ! Slab thickness, m
    dx_tentative    = getParameterValue(25) ! tentative (might differ from definitive) grid size along slab length , m
    dy_tentative    = getParameterValue(26) ! tentative (might differ from definitive) grid size along slab thickness, m
    Slab_half_thickness = Slab_thickness/2.0 ! Slab half thickness, m
    
    M = nint(Slab_length/dx_tentative)              ! Number of divisions along slab length (x direction), -   
    N = nint(Slab_half_thickness/dy_tentative)      ! Number of divisions along half the slab thickness (y direction), - 
    
    Static_storage_size = 4*(M+1) + 4*(M+1)*(N+1)   ! Total static storage spots needed for Ta, Ta_new, Tw, Tw_new, H, H_new, T and T_new
    
    !   Tell the TRNSYS Engine How This Type Works 
    Call SetNumberofParameters(30)            
    Call SetNumberofInputs(3)            
    Call SetNumberofDerivatives(0)            
    Call SetNumberofOutputs(7)            
    Call SetIterationMode(1)     ! 1 indicates the type should be caller every iteration       
    Call SetNumberStoredVariables(Static_storage_size,0) 
    
    !  Set the Correct Input and Output Variable Types 
    Call SetInputUnits(1,'TE1')     ! Inlet air temperature
    Call SetInputUnits(2,'MF2')     ! Mass flow rate
    Call SetInputUnits(3,'DM1')     ! Inlet air humidity ratio
    
    Call SetOutputUnits(1,'TE1')    ! Outlet air temperature 
    Call SetOutputUnits(2,'MF2')    ! Mass flow rate    
    Call SetOutputUnits(3,'DM1')    ! Outlet air humidity ratio
    Call SetOutputUnits(4,'TE1')    ! Average pcm temperature
    Call SetOutputUnits(5,'EN1')    ! Total enthalpy
    Call SetOutputUnits(6,'DM1')    ! Lambda average
    Call SetOutputUnits(7,'DM1')    ! Number of iterations
    Return 
EndIf 


!------------------------------------------------------------------------------ 
    !DO ALL OF THE "START TIME" MANIPULATIONS HERE - There Are No Iterations at the Intial Time 
If (getIsStartTime()) Then
    !  Read in the Values of the Parameters from the Input File 
    relax   = getParameterValue(1) ! Underrelaxation factor (between 0 and 1)
    rho_s   = getParameterValue(2) ! PCM solid density, kg/m3 
    rho_l   = getParameterValue(3) ! PCM liquid density, kg/m3
    ks      = getParameterValue(4) ! PCM solid thermal conductivity, W/m-K
    kl      = getParameterValue(5) ! PCM liquid thermal conductivity, W/m-K
    cs      = getParameterValue(6) ! PCM solid specific heat, J/kg-K
    cl      = getParameterValue(7) ! PCM liquid specific heat, J/kg-K
    L       = getParameterValue(8) ! PCM latent heat, J/kg
    Tm      = getParameterValue(9) ! PCM melting temperature midpoint,  C
    T_range = getParameterValue(10) ! PCM melting temperature range,  C
    rho_a   = getParameterValue(11) ! Air density (considered constant), kg/m3
    c_a     = getParameterValue(12) ! Air specific heat (considered constant), J/kg-K
    mu_a    = getParameterValue(13) ! Air dynamic viscosity (considered constant), Pa-s
    k_a     = getParameterValue(14) ! Air thermal conductivity (considered constant), W/m-K
    Pr      = getParameterValue(15) ! Air Prandtl number
    rho_w   = getParameterValue(16) ! Container wall density, kg/m3
    thickn_w = getParameterValue(17) ! Container wall thickness, m
    c_w     = getParameterValue(18) ! Container wall specific heat, J/kg-K
    k_w     = getParameterValue(19) ! Container wall thermal conductivity, W/m-K
    N_slabs = getParameterValue(20) ! Number of slabs
    Slab_length     = getParameterValue(21) ! Slab length, m
    Slab_width_s      = getParameterValue(22) ! This width refers to the width/height of PCM in slab in solid phase 
    Slab_thickness  = getParameterValue(23) ! Slab thickness, m
    Air_gap         = getParameterValue(24) ! Thickness of air gap between slabs, m
    dx_tentative    = getParameterValue(25) ! tentative (might differ from definitive) grid size along slab length , m
    dy_tentative    = getParameterValue(26) ! tentative (might differ from definitive) grid size along slab thickness, m
    TOL_delta_H     = getParameterValue(27) ! tolerance for error between iterations for PCM node enthalpy
    TOL_delta_Ta    = getParameterValue(28) ! tolerance for error between iterations for air node temperature 
    TOL_delta_Tw    = getParameterValue(29) ! tolerance for error between iterations for wall node temperature
    T0              = getParameterValue(30) ! Initial storage temperature
   
    !   Check the parameters for problems and RETURN if any are found 
    IF((relax.GT.1.).OR.(relax.LE.0.)) CALL foundBadParameter(1,'Fatal','The underrelax factor must be between 0 and 1')
    IF(rho_s.LE.0.)     CALL foundBadParameter(2,'Fatal','The solid density must be positive') 
    IF(rho_l.LE.0.)     CALL foundBadParameter(3,'Fatal','The liquid density must be positive')
    IF(ks.LE.0.)        CALL foundBadParameter(4,'Fatal','The solid conductivity must be positive') 
    IF(kl.LE.0.)        CALL foundBadParameter(5,'Fatal','The liquid conductivity must be positive')
    IF(cs.LE.0.)        CALL foundBadParameter(6,'Fatal','The solid specific heat must be positive')
    IF(cl.LE.0.)        CALL foundBadParameter(7,'Fatal','The liquid specific heat must be positive')
    IF(L.LE.0.)         CALL foundBadParameter(8,'Fatal','The latent heat must be positive') 
    IF(T_range.LE.0.)   CALL foundBadParameter(10,'Fatal','The melting T range must be positive')
    IF(rho_a.LE.0.)     CALL foundBadParameter(11,'Fatal','The air density must be positive')
    IF(c_a.LE.0.)       CALL foundBadParameter(12,'Fatal','The air specific heat must be positive')        
    IF(mu_a.LE.0.)      CALL foundBadParameter(13,'Fatal','The air viscosity must be positive')
    IF(k_a.LE.0.)       CALL foundBadParameter(14,'Fatal','The air conductivity must be positive')
    IF(Pr.LE.0.)        CALL foundBadParameter(15,'Fatal','The air Prandtl number must be positive')
    IF(rho_w.LE.0.)     CALL foundBadParameter(16,'Fatal','The wall density must be positive')
    IF(thickn_w.LE.0.)   CALL foundBadParameter(17,'Fatal','The wall thickness must be positive')
    IF(c_w.LE.0.)       CALL foundBadParameter(18,'Fatal','The wall specific heat must be positive')
    IF(k_w.LE.0.)       CALL foundBadParameter(19,'Fatal','The wall conductivity must be positive')
    IF(N_slabs.LT.1)    CALL foundBadParameter(20,'Fatal','The number of slabs must be a positive integer')
    IF(Slab_length.LE.0.)   CALL foundBadParameter(21,'Fatal','The slab length must be positive')
    IF(Slab_width_s.LE.0.)    CALL foundBadParameter(22,'Fatal','The slab width must be positive')
    IF(Slab_thickness.LE.0.)    CALL foundBadParameter(23,'Fatal','The slab thickness must be positive')
    IF(Air_Gap.LE.0.)           CALL foundBadParameter(24,'Fatal','The air gap thickness must be positive')
    IF(dx_tentative.LE.0.)                CALL foundBadParameter(25,'Fatal','dx must be positive')
    IF(dy_tentative.LE.0.)                CALL foundBadParameter(26,'Fatal','dy must be positive')
    IF(TOL_delta_H.LE.0.)       CALL foundBadParameter(27,'Fatal','TOL_delta_H must be positive')
    IF(TOL_delta_Ta.LE.0.)      CALL foundBadParameter(28,'Fatal','TOL_delta_Ta must be positive')
    IF(TOL_delta_Tw.LE.0.)      CALL foundBadParameter(29,'Fatal','TOL_delta_Tw must be positive')
        
    IF (ErrorFound()) RETURN
    
    M = nint(Slab_length/dx_tentative)              ! Number of divisions along slab length (x direction), -   
    N = nint(Slab_half_thickness/dy_tentative)      ! Number of divisions along half the slab thickness (y direction), - 
    
    !SET INITIAL VALUES OF STORAGE VARIABLES
    
    ! Calculating value of H0 based on T0 
    Tl = Tm + T_range/2.         !Higher end of melting T range
    Ts = Tm - T_range/2.         !Lower end of melting T range
        
    if (T0 <= Ts) then           
        H0 = cs*(T0 - Ts)
    else if (T0 >= Tl) then
        H0 = L + cl*(T0 - Tl)
    else
        H0 = L*((T0 - Ts)/T_range)
    end if           
    
    !setting initial values for Ta in storage array
    i = 0
    do j = 1,M+1
        call setStaticArrayValue(i*(M+1) + j,T0)
    end do
        
    !setting initial values for Tw in storage array
    i = 1
    do j = 1,M+1
        call setStaticArrayValue(i*(M+1) + j,T0)
    end do
        
    !setting initial values for T in storage array
    do i = 2, N+2
        do j = 1,M+1
            call setStaticArrayValue(i*(M+1) + j,T0)
        end do
    end do
              
    !setting initial values for H in storage array
    do i = N+3, 2*N+3
        do j = 1,M+1
            call setStaticArrayValue(i*(M+1) + j,H0)
        end do
    end do
        
    !setting initial values for Ta_new in storage array
    i = 2*N+4
    do j = 1,M+1
        call setStaticArrayValue(i*(M+1) + j,T0)
    end do

    !setting initial values for Tw_new in storage array
    i = 2*N+5
    do j = 1,M+1
        call setStaticArrayValue(i*(M+1) + j,T0)
    end do
    
    !setting initial values for T_new in storage array
    do i = 2*N+6, 3*N+6
        do j = 1,M+1
            call setStaticArrayValue(i*(M+1) + j,T0)
        end do
    end do
    
    !setting initial values for H_new in storage array
    do i = 3*N+7, 4*N+7
        do j = 1,M+1
            call setStaticArrayValue(i*(M+1) + j,H0)
        end do
    end do
    
    grid_weights = 0            ! First set all elements to 0
    do i = 0,N                            ! Now set all elements in the actual matrix to 1
       do j = 0,M
            grid_weights(i,j) = 1.0
        end do
    end do
        
    j=0                         ! Here we assign the corresponding weights to the control volumes.  
    do i = 0,N-1
        grid_weights(i,j) = 0.5   
    end do
    
    j=M    
    do i = 0,N-1
        grid_weights(i,j) = 0.5   
    end do
    
    grid_weights(N,0) = 0.25
    grid_weights(N,M) = 0.25
    
    i=N
    do j = 1,M-1
        grid_weights(i,j) = 0.5
    end do
    
    if (H0 <= Hs) then
        Lambda = 0
    else if (H0 >= Hl) then
        Lambda = 1.0
    else
        Lambda = H0/Hl
    end if
    
    Mass_pcm = N_slabs*(Slab_length*Slab_width_s*Slab_thickness)*rho_s     ! Total mass of PCM in LHS, kg
    Lambda_avg = sum(Lambda*grid_weights)/sum(grid_weights) ! Calculates the weighted average liquid fraction of the entire PCM region
    Total_enthalpy = H0*Mass_pcm/1000   ! Calculates the total enthalpy stored from the lower end of the phase change, kJ
       
    !  Set the Initial Values of the Outputs
    Call setOutputValue(1, T0)
    Call setOutputValue(2, getInputValue(2))
    Call setOutputValue(3, getInputValue(3))
    Call setOutputValue(4, T0)
    Call setOutputValue(5, Total_enthalpy)
    Call setOutputValue(6, Lambda_avg)
    Call setOutputValue(7, 0)
    
    RETURN 
EndIf
     
!------------------------------------------------------------------------------ 
! "MULTIPLE UNIT" Manipulations ReRead the Parameters if Another Unit of This Type Has Been Called Last 
If (getIsReReadParameters()) Then 
    relax   = getParameterValue(1) ! Underrelaxation factor (between 0 and 1)
    rho_s   = getParameterValue(2) ! PCM solid density, kg/m3 
    rho_l   = getParameterValue(3) ! PCM liquid density, kg/m3
    ks      = getParameterValue(4) ! PCM solid thermal conductivity, W/m-K
    kl      = getParameterValue(5) ! PCM liquid thermal conductivity, W/m-K
    cs      = getParameterValue(6) ! PCM solid specific heat, J/kg-K
    cl      = getParameterValue(7) ! PCM liquid specific heat, J/kg-K
    L       = getParameterValue(8) ! PCM latent heat, J/kg
    Tm      = getParameterValue(9) ! PCM melting temperature midpoint,  C
    T_range = getParameterValue(10) ! PCM melting temperature range,  C
    rho_a   = getParameterValue(11) ! Air density (considered constant), kg/m3
    c_a     = getParameterValue(12) ! Air specific heat (considered constant), J/kg-K
    mu_a    = getParameterValue(13) ! Air dynamic viscosity (considered constant), Pa-s
    k_a     = getParameterValue(14) ! Air thermal conductivity (considered constant), W/m-K
    Pr      = getParameterValue(15) ! Air Prandtl number
    rho_w   = getParameterValue(16) ! Container wall density, kg/m3
    thickn_w = getParameterValue(17) ! Container wall thickness, m
    c_w     = getParameterValue(18) ! Container wall specific heat, J/kg-K
    k_w     = getParameterValue(19) ! Container wall thermal conductivity, W/m-K
    N_Slabs = getParameterValue(20) ! Number of slabs
    Slab_length     = getParameterValue(21) ! Slab length, m
    Slab_width_s      = getParameterValue(22) ! This width refers to the width/height of PCM in slab in solid phase 
    Slab_thickness  = getParameterValue(23) ! Slab thickness, m
    Air_gap         = getParameterValue(24) ! Thickness of air gap between slabs, m
    dx_tentative    = getParameterValue(25) ! tentative (might differ from definitive) grid size along slab length , m
    dy_tentative    = getParameterValue(26) ! tentative (might differ from definitive) grid size along slab thickness, m
    TOL_delta_H     = getParameterValue(27) ! tolerance for error between iterations for PCM node enthalpy
    TOL_delta_Ta    = getParameterValue(28) ! tolerance for error between iterations for air node temperature 
    TOL_delta_Tw    = getParameterValue(29) ! tolerance for error between iterations for wall node temperature
    T0              = getParameterValue(30) ! Initial storage temperature
EndIf

!----------------------------------------------------------------------
! "EVERY TIME STEP" manipulations 

! RETRIEVE STORED VALUES

! Retrieve Ta array
i=0 
do j = 1,M+1 
    Ta(j-1) = getStaticArrayValue(i*(M+1) + j)
end do    

! Retrieve Tw array    
i=1
do j = 1,M+1    
    Tw(j-1) = getStaticArrayValue(i*(M+1) + j)
end do

! Retrieve T array
do i=2,N+2
    do j=1,M+1
        T(i-2,j-1) = getStaticArrayValue(i*(M+1) + j)    
    end do    
end do

! Retrieve H array
do i=N+3,2*N+3
    do j=1,M+1
        H(i-(N+3),j-1) = getStaticArrayValue(i*(M+1) + j)    
    end do    
end do

! Retrieve Ta_new array
i=2*N+4
do j = 1,M+1 
    Ta_new(j-1) = getStaticArrayValue(i*(M+1) + j)
end do

! Retrieve Tw_new array
i=2*N+5
do j = 1,M+1 
    Tw_new(j-1) = getStaticArrayValue(i*(M+1) + j)
end do
 
! Retrieve T_new array
do i=2*N+6,3*N+6
    do j=1,M+1
        T_new(i-(2*N+6),j-1) = getStaticArrayValue(i*(M+1) + j)    
    end do    
end do

! Retrieve H_new array
do i=3*N+7,4*N+7
    do j=1,M+1
        H_new(i-(3*N+7),j-1) = getStaticArrayValue(i*(M+1) + j)    
    end do    
end do

!Retrieve Current Inputs to the Model 
Tair_in = getInputValue(1)    
mdot_a = getInputValue(2)
Xa_in = getInputValue(3)
    
If (mdot_a <= 0) Call foundBadInput(2,'Fatal', 'The input flow rate must be positive')    
If (ErrorFound() ) Return

! Perform calculations
    
! First constant values derived from parameters are calculated
expansion = rho_s/rho_l    ! Volume expansion from completely solid to completely melted
Hs = 0.0                    !Enthalpy of solid at Tm 
Hl = Hs + L                !Enthalpy of liquid at Tm

Tl = Tm + T_range/2.         !Higher end of melting T range
Ts = Tm - T_range/2.         !Lower end of melting T range
N_channels = N_slabs    ! Assuming that there is a channel between inner slabs and the outer channels adjacent to the lateral walls are half as wide.
Duct_width = Slab_width_s*expansion  ! This would be equal to the slab width when the whole PCM is liquid, due to the expansion, &
                                       ! and would be equal to the duct width, m
Slab_half_thickness = Slab_thickness/2.0  ! The half thickness, which is the actual modeled region, m
   
Mass_pcm = N_slabs*(Slab_length*Slab_width_s*Slab_thickness)*rho_s     ! Total mass of PCM in LHS, kg

Half_air_gap = Air_gap/2.                   ! Half thickness, which is what is taken in the model since we are considering only the ehat transfer to/from one slab and using symmetry

Dh = 4*Duct_width*Air_gap/(2*(Duct_width + Air_gap))    ! Hydraulic diameter of air channels, assuming the slabs go
                                                            ! all the way to the walls
                                                            
!Gap_AR = Air_gap/Duct_width                            ! Aspect ratio of the air gaps
    
A = 1/(rho_a*c_a*Half_air_gap)        ! Factor used in equation for the air to simplify writing
    

M = nint(Slab_length/dx_tentative)              ! Number of divisions along slab length (x direction), -   
N = nint(Slab_half_thickness/dy_tentative)      ! Number of divisions along half the slab thickness (y direction), - 
dy = Slab_half_thickness/(N+0.5)
    
x_position = 0        ! Elements of vector of node coords. along the slab length, m, initialized to 0 
y_position = 0        ! Elements of vector of node coords. along the slab half thickness, m, initialized to 0

do j=0,M                        ! Assign values of x_position according to slab length and M
    x_position(j) = j*Slab_length/M
end do

do i=0,N                        ! Assign values of y_position according to slab half thickness and N
    y_position(i) = dy/2+i*dy
end do
   
dx = x_position(1) - x_position(0)         ! Mesh size along slab lenth (actual), m
dy = y_position(1) - y_position(0)         ! Mesh size along slab thickness (actual), m

fe = (dy/2.)/(dy/2. + thickn_w)         ! ratio of distances for harmonic mean interface conductivity between container 
                                            ! wall and PCM (Libro Patankar page. 45)
!Store_dt = Store_dt_min*60                  ! Data storage interval, s
!Store_dt_details = Store_dt_min_details*60  ! Data storage interval for detailed slab data, s
    
! Now come the definition of variables and arrays that vary with iterations and timesteps
V_a = mdot_a/rho_a
va =  V_a/(N_channels*Air_gap*Duct_width)              ! Mean air velocity through gaps, m/s
Re = rho_a*va*Dh/mu_a
ff = (0.79*log(Re) - 1.64)**(-2)   ! Friction factor to be used in Nu calculation with Gnielinski eq.
Xaster = (Slab_length/Dh)/(Re*Pr)
if (Re <= 2300.0) then
    Nu = 7.55 + (0.024*Xaster**(-1.14))/(1 + 0.0358*(Pr**0.17)*Xaster**(-0.64)) ! Nu EN PARALLEL PLATES SIMULTANEAOUSLY DEVELOPING FLOW ROHSENOW PAG 5.63
    !Nu = 7.541*(1 - 2.61*Gap_AR + 4.97*Gap_AR**2 - 5.119*Gap_AR**3 + 2.702*Gap_AR**4 - 0.548*Gap_AR**5)
else
    Nu = (ff/8.)*(Re-1000.)*Pr/(1 + 12.7*(ff/8.)**0.5*(Pr**(2./3.)-1.))    ! Nusselt number by Gnielinski (appears in &
                                                                               ! practically all sources as a good eq., used by Dolado and others)
end if

h_aw = Nu*k_a/Dh

delta_Ta = 0        !Matrix of delta_Ta initialized to 0 (from j=0 to M_max)
delta_Tw = 0         !Matrix of delta_Tw initialized to 0(from j=0 to M_max)
delta_H = 0         !Matrix of delta_H initialized to 0 (from i=0 to N_max and j=0 to M_max)

do j = 0,M              !Set to 1 Values of delta_Ta, delta_Tw and delta_H where where our geometry actually exists (i=0 to N and j=0 to M)
    delta_Ta(j) = 1.0
    delta_Tw(j) = 1.0
end do

do i = 0,N
    do j = 0,M
        delta_H(i,j) = 1.0 
    end do
end do


do i = 0,N              ! Setting values of matrices k and Lambda according to H
    do j = 0,M
        if (H(i,j) <= Hs) then
            k(i,j) = ks
            Lambda(i,j) = 0
        else if (H(i,j) >= Hl) then
            k(i,j) = kl
            Lambda(i,j) = 1.0
        else
            k(i,j) = ks + (H(i,j) - Hs)/(Hl - Hs)*(kl - ks) 
            Lambda(i,j) = H(i,j)/Hl
        end if
    end do
end do

grid_weights = 0            ! First set all elements to 0
do i = 0,N                            ! Now set all elements in the actual matrix to 1
    do j = 0,M
        grid_weights(i,j) = 1.0
    end do
end do
        
j=0                         ! Here we assign the corresponding weights to the control volumes.  
do i = 0,N-1
    grid_weights(i,j) = 0.5   
end do
    
j=M    
do i = 0,N-1
    grid_weights(i,j) = 0.5   
end do
    
grid_weights(N,0) = 0.25
grid_weights(N,M) = 0.25
    
i=N
do j = 1,M-1
    grid_weights(i,j) = 0.5
end do

iterations = 0
dt = dt_hours*3600    
        
do while (maxval(delta_H) > TOL_delta_H .or. maxval(delta_Ta) > TOL_delta_Ta .or. maxval(delta_Tw) > TOL_delta_Tw)
    ! LOOP OVER AIR NODES FROM i=0 TO M 
    Ta_iter = Tair_in
    Ta_new(0) = Ta_iter
    delta_Ta(0) = Ta_iter - Ta_new(0)
            
    do j = 1,M
        Ta_iter = (Ta(j) + (va*dt/dx)*Ta_new(j-1) + dt*A*h_aw*Tw_new(j))/(1 + dt*(va/dx + A*h_aw))
        delta_Ta(j) = abs(Ta_iter - Ta_new(j))
        Ta_new(j) = Ta_iter
    end do
            
    ! LOOP OVER CONTAINER WALL NODES
            
    do j = 0,M
                
        if (j>0) then ! This is to be able to use the same equations for internal and boundary node
            jm1 = j-1 
        else
            jm1 = j+1
        end if
                
        if (j<M) then 
            jp1 = j+1 
        else 
            jp1 = j-1
        end if
                
        k_int = 1/((1-fe)/k_w + fe/(k(0,j)))  ! Harmonic mean interface conducticvity for interface between container wall &
                                                ! and PCM (Patankar). k[0,j] is the conduct of the PCM adjacent to wall
                                                      
        W = rho_w*thickn_w*c_w + (2*dt*k_w*thickn_w)/(dx**2) + dt*h_aw + dt*k_int/(thickn_w + dy/2) ! Denominator in equation &
                                                                                                        !for container wall
                                                                                                            
        Tw_iter = (rho_w*thickn_w*c_w*Tw(j) + (dt*k_w*thickn_w/dx**2)*(Tw_new(jp1) + Tw_new(jm1)) + dt*h_aw*Ta_new(j) + &
            (dt*k_int/(thickn_w + dy/2))*T_new(0,j))/W

        Tw_iter = relax*Tw_iter + (1-relax)*Tw_new(j)             !Applying underrelaxation 

        delta_Tw(j) = abs(Tw_iter - Tw_new(j))
                
        Tw_new(j) = Tw_iter
    end do
            
    ! LOOP OVER THE NODES OF THE SURFACE ROW (i=0) FROM [0,0] UNTIL [0,M]
    i = 0
    do j = 0,M
        if (j>0) then ! This is to be able to use the same equations for internal and boundary node
            jm1 = j-1 
        else
            jm1 = j+1
        end if
                
        if (j<M) then 
            jp1 = j+1 
        else 
            jp1 = j-1 
        end if
               
        k_int = 1/((1-fe)/k_w + fe/(k(0,j)))
                
        kx_ant = 2*k(i,jm1)*k(i,j)/(k(i,jm1)+k(i,j))  !Patankar formulas for interface conductivity using &
        kx_post = 2*k(i,jp1)*k(i,j)/(k(i,jp1)+k(i,j))  !harmonic mean
        ky_post = 2*k(i+1,j)*k(i,j)/(k(i+1,j)+k(i,j))      
               
        if (H_new(i,j) <= Hs) then
            S = rho_s/dt
        else if (Hs < H_new(i,j) .and. H_new(i,j) < Hl) then
            S = ((rho_l - rho_s)/(Hl - Hs)*H_new(i,j) + rho_s)/dt
        else
            S = rho_l/dt
        end if
                
        R = S*H(i,j) + k_int/((thickn_w + dy/2)*dy)*Tw_new(j) + (kx_post/dx**2)*T_new(i,jp1) + &
            (kx_ant/dx**2)*T_new(i,jm1) + (ky_post/dy**2)*T_new(i+1,j)
        F = (k_int/((thickn_w + dy/2)*dy) + kx_post/dx**2 + kx_ant/dx**2 + ky_post/dy**2)
                
        H_iter = (R + F*(Hs/cs - Ts))/(S + F/cs)
        H_iter = relax*H_iter + (1-relax)*H_new(i,j)             ! Applying underrelaxation
        T_iter = (H_iter - Hs)/cs + Ts

        if (Hs <= H_iter .and. H_iter <= Hl) then
            H_iter = (R + F*Hs*T_range/L - F*Ts)/(S + F*T_range/L)
            H_iter = relax*H_iter + (1-relax)*H_new(i,j)            ! Applying underrelaxation
            T_iter = Ts + (H_iter - Hs)*T_range/L
        end if
                
        if (H_iter > Hl) then
            H_iter = (R + F*(Hl/cl - Tl))/(S + F/cl)
            H_iter = relax*H_iter + (1-relax)*H_new(i,j)             ! Applying underrelaxation
            T_iter = (H_iter - Hl)/cl + Tl
        end if
                
        delta_H(i,j) = abs(H_iter - H_new(i,j))  ! Recording the change in enthalpy from past iteration
            
                
        H_new(i,j) = H_iter ! Updating value of H_new[i,j] and T_new[i,j]
        T_new(i,j) = T_iter
            
        if (H_new(i,j) <= Hs) then
            k(i,j) = ks
        else if (Hs < H_new(i,j) .and. H_new(i,j) < Hl) then
            k(i,j) = ks + (kl - ks)/(Hl - Hs)*(H_new(i,j) - Hs) 
        else
            k(i,j) = kl
        end if
    end do
            
    !NOW WE MOVE TO THE INNER ROWS, FROM ROW i=1 TO i=N-1, AND WITHIN EACH ROW MOVING FROM COLUMN j=0 TO j=M
    do i = 1,N-1
        do j = 0,M
            if (j>0) then ! This is to be able to use the same equations for internal and boundary node
                jm1 = j-1 
            else
                jm1 = j+1
            end if
              
            if (j<M) then 
                jp1 = j+1 
            else 
                jp1 = j-1 
            end if 
                    
            kx_ant = 2*k(i,jm1)*k(i,j)/(k(i,jm1)+k(i,j))  !Patankar formulas for interface conductivity using &
            kx_post = 2*k(i,jp1)*k(i,j)/(k(i,jp1)+k(i,j)) !harmonic mean
            ky_ant = 2*k(i-1,j)*k(i,j)/(k(i-1,j)+k(i,j))  
            ky_post = 2*k(i+1,j)*k(i,j)/(k(i+1,j)+k(i,j)) 
            
            if (H_new(i,j) <= Hs) then
                S = rho_s/dt
            else if (Hs < H_new(i,j) .and. H_new(i,j) < Hl) then
                S = ((rho_l - rho_s)/(Hl - Hs)*H_new(i,j) + rho_s)/dt
            else
                S = rho_l/dt
            end if
                    
            R = S*H(i,j) + (kx_post/dx**2)*T_new(i,jp1) + (kx_ant/dx**2)*T_new(i,jm1) + &
                (ky_post/dy**2)*T_new(i+1,j) + (ky_ant/dy**2)*T_new(i-1,j)
                
            F = (kx_post/dx**2 + kx_ant/dx**2 + ky_post/dy**2 + ky_ant/dy**2)
                    
            H_iter = (R + F*(Hs/cs - Ts))/(S + F/cs)
            H_iter = relax*H_iter + (1-relax)*H_new(i,j)             ! Applying underrelaxation
            T_iter = (H_iter - Hs)/cs + Ts

            if (Hs <= H_iter .and. H_iter <= Hl) then
                H_iter = (R + F*Hs*T_range/L - F*Ts)/(S + F*T_range/L) 
                H_iter = relax*H_iter + (1-relax)*H_new(i,j)             ! Applying underrelaxation
                T_iter = Ts + (H_iter - Hs)*T_range/L
            end if
                    
            if (H_iter > Hl) then
                H_iter = (R + F*(Hl/cl - Tl))/(S + F/cl)
                H_iter = relax*H_iter + (1-relax)*H_new(i,j)             ! Applying underrelaxation
                T_iter = (H_iter - Hl)/cl + Tl
            end if

            delta_H(i,j) = abs(H_iter - H_new(i,j))    ! Recording the change in enthalpy from past iteration

                    
            H_new(i,j) = H_iter    ! Updating value of H_new[i,j] and T_new[i,j]
            T_new(i,j) = T_iter
                                    
            if (H_new(i,j) <= Hs) then
                k(i,j) = ks     
            else if (Hs < H_new(i,j) .and. H_new(i,j) < Hl) then
                k(i,j) = ks + (kl - ks)/(Hl - Hs)*(H_new(i,j) - Hs)
            else
                k(i,j) = kl
            end if
        end do
    end do
            
    ! NOW THE SYMMETRY ROW i=N, GOING FROM COLUMN J=0 TO J
    i = N
    do j = 0,M
        if (j>0) then ! This is to be able to use the same equations for internal and boundary node
            jm1 = j-1 
        else
            jm1 = j+1
        end if
                
        if (j<M) then 
            jp1 = j+1 
        else 
            jp1 = j-1 
        end if
                
        kx_ant = 2*k(i,jm1)*k(i,j)/(k(i,jm1)+k(i,j))   ! Patankar formulas for interface conductivity using &
        kx_post = 2*k(i,jp1)*k(i,j)/(k(i,jp1)+k(i,j))  ! harmonic mean
        ky_ant = 2*k(i-1,j)*k(i,j)/(k(i-1,j)+k(i,j))
                
        if (H_new(i,j) <= Hs) then
            S = rho_s/dt
        else if (Hs < H_new(i,j) .and. H_new(i,j) < Hl) then
            S = ((rho_l - rho_s)/(Hl - Hs)*H_new(i,j) + rho_s)/dt
        else
            S = rho_l/dt
        end if
                
        R = S*H(i,j) + (kx_post/dx**2)*T_new(i,jp1) + (kx_ant/dx**2)*T_new(i,jm1) + &
            2*(ky_ant/dy**2)*T_new(i-1,j)
    
        F = (kx_post/dx**2 + kx_ant/dx**2 + 2*ky_ant/dy**2)
                
        H_iter = (R + F*(Hs/cs - Ts))/(S + F/cs)
        H_iter = relax*H_iter + (1-relax)*H_new(i,j)     ! Applying underrelaxation
        T_iter = (H_iter - Hs)/cs + Ts 
                
        if (Hs <= H_iter .and. H_iter <= Hl) then
            H_iter = (R + F*Hs*T_range/L - F*Ts)/(S + F*T_range/L) 
            H_iter = relax*H_iter + (1-relax)*H_new(i,j)             ! Applying underrelaxation
            T_iter = Ts + (H_iter - Hs)*T_range/L
        end if
                
        if (H_iter > Hl) then
            H_iter = (R + F*(Hl/cl - Tl))/(S + F/cl)
            H_iter = relax*H_iter + (1-relax)*H_new(i,j)             ! Applying underrelaxation
            T_iter = (H_iter - Hl)/cl + Tl     
        end if
                
        delta_H(i,j) = abs(H_iter - H_new(i,j))
                
        H_new(i,j) = H_iter
        T_new(i,j) = T_iter
               
        if (H_new(i,j) <= Hs) then
            k(i,j) = ks
        else if (Hs < H_new(i,j) .and. H_new(i,j) < Hl) then
            k(i,j) = ks + (kl - ks)/(Hl - Hs)*(H_new(i,j) - Hs)
        else
            k(i,j) = kl
        end if
    end do
            
    iterations = iterations + 1
    !iterations = 34    
        
end do

do i = 0,N              ! Calculating the values of matrix Lamdba at the new time
    do j = 0,M
        if (H_new(i,j) <= Hs) then
            Lambda(i,j) = 0
        else if (H_new(i,j) >= Hl) then
            Lambda(i,j) = 1.0
        else
            Lambda(i,j) = H_new(i,j)/Hl
        end if
    end do
end do

            
pcm_width = Slab_width_s + Lambda_avg*(Duct_width - Slab_width)   ! This PCM width is the actual width/height of the pcm inside
                                                ! the slabs as it melts according to how much pcm has melted (liquid fraction).
                                                ! it can go from the level when all pcm is solid, to the level when all pcm is 
                                                ! liquid, according to the calculated average liquid fraction
        
Ta_out_mixed = (Ta_new(M)*pcm_width + Tair_in*(Duct_width - pcm_width))/Duct_width    ! This calculates the mixed
                                                ! temperature of the air at the outlet of the slab and the air which was
                                                ! above the pcm level,which leaves de LHS as it arrived, with the inlet T
                                                ! which is the collector outlet temperature
        
!Tout_avg_stor_dt =  (Tout_avg_stor_dt*(Count_stor_dt - 1) + Ta_out_mixed)/Count_stor_dt    ! This is used to calculate the 
                                                ! average outlet T during the store_dt. It calculates the mean of all Tout values
                                                ! between store_dts to give a value of mean Tout for the whole storage period. 
                                                ! For example if i store data every minute and my ts is 10s, the code calculates 
                                                ! a Tout every 10s, so here we calculate the mean of 6 Tout values which represent 
                                                ! the mean Tout of that minute. 
   
! SET STORAGE VALUES

!setting provisional end of timestep values for Ta_new in storage array
i = 2*N+4
do j = 1,M+1
    call setStaticArrayValue(i*(M+1) + j,Ta_new(j-1))
end do

!setting provisional end of timestep values for Tw_new in storage array
i = 2*N+5
do j = 1,M+1
    call setStaticArrayValue(i*(M+1) + j,Tw_new(j-1))
end do
    
!setting provisional end of timestep values for T_new in storage array
do i = 2*N+6, 3*N+6
    do j = 1,M+1
        call setStaticArrayValue(i*(M+1) + j,T_new(i-(2*N+6),j-1))
    end do
end do
    
!setting provisional end of timestep values for H_new in storage array
do i = 3*N+7, 4*N+7
    do j = 1,M+1
       call setStaticArrayValue(i*(M+1) + j,H_new(i-(3*N+7),j-1))
    end do
end do

Lambda_avg = sum(Lambda*grid_weights)/sum(grid_weights) ! Calculates the weighted average liquid fraction of the entire PCM region
Total_enthalpy = sum(H_new*grid_weights)/sum(grid_weights)*Mass_pcm/1000   ! Calculates the total enthalpy stored from the lower end of the phase change
Avg_temperature = sum(T_new*grid_weights)/sum(grid_weights)  ! Calculates the weighted average PCM temperature

! SET OUTPUT VALUES    
Call setOutputValue(1, Ta_out_mixed)
Call setOutputValue(2, mdot_a)
Call setOutputValue(3, Xa_in)
Call setOutputValue(4, Avg_temperature)
Call setOutputValue(5, Total_enthalpy)
Call setOutputValue(6, Lambda_avg)  
Call setOutputValue(7, iterations)
    

Return
End Subroutine TYPE201