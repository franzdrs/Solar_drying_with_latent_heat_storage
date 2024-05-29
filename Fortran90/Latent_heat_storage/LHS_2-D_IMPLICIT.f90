module file_data_module
    implicit none
    integer :: num_rows
    character(len=12), allocatable :: Date(:), Hour_min(:)
    double precision, allocatable :: Hour_dec(:), t_min(:), Ins(:), Tamb(:), wind(:), Patm(:), X(:), Ta_in(:), mdot_a(:)
end module file_data_module

! SUBROUTINE TO READ THE INPUT DATA FILE AND ASSIGN EACH DATA VARIABLE TO A LIST WHOSE LENGTH IS THE NUMBER OF DATA ROWS.
! AS IT IS, THE SUBROUTINE IS MADE TO READ THE OUPUT FILE FROM THE SAH PROGRAM WHICH HAS 6 ROWS WITH SETUP INFO, WHICH ARE SKIPPED
! THE USER CAN MODIFY THE SUBROUTINE OR THE INPUT FILE TO MATCH IF THE DATA FILE FORMAT IS DIFFERENT
subroutine read_file(input_file)
    use file_data_module
    implicit none
    character(len=*), intent(in) :: input_file
    integer :: input_unit_number, i, ierr
    character(len=150) :: line

    input_unit_number = 10
    
    !Open the file
    open(unit=input_unit_number, file=input_file, status='old')

    !Skip the first 6 rows
    do i = 0, 5
        read(input_unit_number, *)
    end do
    
    !Initialize variables
    num_rows = 0
    ierr = 0

    !First pass: Count the number of rows
    do
        read(input_unit_number, '(A)', iostat=ierr)
        if (ierr /= 0) exit  ! Exit loop if end-of-file or error
        num_rows = num_rows + 1
    end do

    !Allocate memory for the arrays based on the number of rows
    allocate(Date(0:num_rows-1))
    allocate(Hour_min(0:num_rows-1))
    allocate(Hour_dec(0:num_rows-1))
    allocate(t_min(0:num_rows-1))
    allocate(Ins(0:num_rows-1))
    allocate(Tamb(0:num_rows-1))
    allocate(wind(0:num_rows-1))
    allocate(Patm(0:num_rows-1))
    allocate(X(0:num_rows-1))
    allocate(Ta_in(0:num_rows-1))
    allocate(mdot_a(0:num_rows-1))
    
    !Rewind the file to the beginning
    rewind(input_unit_number)

    !Skip the first 6 rows again
    do i = 0, 5
        read(input_unit_number, *)
    end do
    
    !Second pass: Read data and store in arrays
    do i = 0, num_rows-1
        read(input_unit_number, '(A)') line
        !print *, line
        !Parse the line and store values in arrays
        read(line, *) Date(i), Hour_min(i), Hour_dec(i), t_min(i), Ins(i), Tamb(i), wind(i), Patm(i), X(i), Ta_in(i), mdot_a(i)
    end do

    !Close the file
    close(input_unit_number)
    
end subroutine read_file

!HERE STARTS THE MAIN PROGRAM TO CALCULATE THE OPERATION OF THE LATENT HEAT STORAGE    
program LHS_2D
    use file_data_module
    implicit none

    double precision :: relax
    double precision :: rho_s, rho_l, expansion, Mass_pcm, ks, kl, cs, cl, Hs, L, Hl, Ts, T_range, Tl, Tm, Lambda_avg, Tair_in
    double precision :: rho_a, c_a, mu_a, k_a, V_a, Pr, mdot_air
    double precision :: rho_w, thickn_w, c_w, k_w, k_int
    double precision :: Slab_length, Slab_width, Slab_thickness, Slab_half_thickness, Air_gap, Half_air_gap, &
        Dh, Duct_width, pcm_width
    double precision :: va, Re, Nu, ff, Tw0, T0, Ta0, h_aw, A, dt, dx, dy, fe, W, R, F, Xaster
    integer :: N_slabs, N_channels
    integer :: M, N, i, j, jm1, jp1
    double precision, allocatable :: x_position(:), y_position(:)
    double precision, allocatable :: Ta(:), Ta_new(:) 
    double precision, allocatable :: Tw(:), Tw_new(:)
    double precision, allocatable :: T(:,:), T_new(:,:)
    double precision, allocatable :: H(:,:), H_new(:,:), Lambda(:,:)
    double precision :: Ta_iter, Tw_iter, T_iter, H_iter, kx_ant, kx_post, ky_ant, ky_post, S
    double precision :: Store_dt_min, Store_dt_min_details, Store_dt, Store_dt_details, Total_sim_t, Curr_sim_time, Data_time_step
    double precision :: Next_store_time, Next_store_time_details, Count_stor_dt, Tout_avg_stor_dt, Ta_out_mixed
    double precision, allocatable :: c(:,:), k(:,:)

    double precision, allocatable :: delta_H(:,:), delta_Ta(:), delta_Tw(:) 
    double precision :: TOL_delta_H, TOL_delta_Ta, TOL_delta_Tw
    integer :: indx, iterations, iterations_tot, output_unit_number,output_DETAILS_T_unit_number, output_DETAILS_H_unit_number
   
    character(len=100) :: output_file = 'LHS_output_Stutt_Aug2022.txt'
    character(len=100) :: output_DETAILS_T_file = 'LHS_output_Stutt_Aug2022_DETAILS_T.txt'
    character(len=100) :: output_DETAILS_H_file = 'LHS_output_Stutt_Aug2022_DETAILS_H.txt'
    
    output_unit_number = 11
    output_DETAILS_T_unit_number = 12
    output_DETAILS_H_unit_number = 13
    
    Data_time_step = 60         ! Time resolution of the data in the input file (s)
    relax = 1.0                 ! Underrelaxation factor, can be lower than 1 for stability if required, but more than 0
    
    !Call the subroutine to read data
    call read_file("Output_SAH_Stutt_Aug2022_as_INPUT.txt")
    
    ! USER INPUTS    
    
    !PCM PROPERTIES
    rho_s = 1517.               ! solid PCM density (kg/m3)
    rho_l = 1442.               ! liquid PCM density (kg/m3)
     
    ks = 1.09                   ! solid PCM conductivity (W/m-K)
    kl = 0.54                   ! liquid PCM conductivity (W/m-K)
                 
    cs = 2000.                  ! solid PCM specific heat, (J/kg-K)
    cl = 2000.                  ! liquid PCM specific heat (J/kg-K)

    Tm = 40.0                   ! Midpoint of PCM melting temperature range ( C)
    T_range = 4.0               ! Temperature range over which the PCM melts ( C)
    
    Hs = 0.0                    ! Enthalpy of solid at Ts (J/kg) (we assume that PCM enthalpy is 0 at Ts so negative below Ts)

    L = 190000.0                ! PCM latent heat (J/kg)    
    Hl = Hs + L                 ! Enthalpy of liquid at Tl (J/kg)
    
    !AIR PROPERTIES ASSUMED CONSTANT, THE SAME AS IN SOLAR AIR HEATER MODEL
    rho_a = 1.127               ! Air density (at about 40  C) (kg/m3)
    c_a = 1007.                 ! Air specific heat (J/kg-K)
    mu_a = 1.907e-5             ! Air viscosity (at 40  C) (Pa-s)
    k_a = 0.02735               ! Themal conductivity of air (at 40  C) (W/m-K)
    Pr = 0.706                  ! Prandtl number of air

    !PROPERTIES OF CONTAINER WALLS
    rho_w = 2670.               ! Container wall density (kg/m3)
    thickn_w = 0.0005           ! Container wall thickness (m)
    c_w = 890.                  ! Container wall specific heat (J/kg-K)
    k_w = 130.                  ! Container wall thermal conductivity (W/m-K)
    
    !STORAGE VESSEL PROPERTIES
    N_slabs = 29                        ! Number of slabs in storage vessel
    N_channels = N_slabs                ! Assuming that there is a channel between inner slabs and the outer channels adjacent to the lateral walls are half as wide. 
    Slab_length = 1.2                   ! Length of PCM slab, m
    Slab_width = 0.5                    ! Height of PCM in the slab when all PCM is solid (m)
    
    Slab_thickness = 0.01               ! Thickness of the slab (m)
    Air_gap = 0.01                      ! Thickness of the air channels between slabs (m)
    
    ! SPACE AND TIME GRID SIZE
    dt = 30.                       ! Timestep (s)
    dx = 0.02                      ! Mesh size along slab length (tentative) (m)
    dy = 0.001                     ! Mesh size along slab thickness (tentative) (m)
    
    ! TOLERANCE FOR THE DIFFERENCE BETWEEN CONSECUTIVE ITERATIONS
    TOL_delta_H = 0.1               ! Tolerance of nodal enthalpy (J/kg)
    TOL_delta_Ta = 0.0001           ! Tolerance of nodal air temperature ( C)
    TOL_delta_Tw = 0.0001           ! Tolerance of nodal wall temperature ( C)
    
    ! INITIAL CONDITIONS
    Tw0 = Tamb(0)                  ! Initial temperature of container wall ( C)
    T0 = Tamb(0)                   ! Initial temperature of PCM ( C)
    Ta0 = Tamb(0)                  ! Initial temperature of air inside storage ( C)

    ! CHOOSE HOW OFTEN TO STORE DATA
    Store_dt_min = 1.    ! Data storage interval (min)
    Store_dt_min_details = 5.    ! Data storage interval for detailed slab data (min)

    ! CALCULATIONS START HERE
    expansion = rho_s/rho_l     
    Tl = Tm + T_range/2         ! Upper bound of melting temperature range ( C)
    Ts = Tm - T_range/2         ! Lower bound of melting temperature range ( C)
    
    Duct_width = Slab_width*expansion   ! This would be equal to the PCM height in the slab when the whole PCM is liquid
                                        ! due to the expansion, and would be equal to the duct width (m)
    Slab_half_thickness = Slab_thickness/2.0  ! The half thickness, which is the actual modeled region (m)
    
    Mass_pcm = N_slabs*(Slab_length*Slab_width*Slab_thickness)*rho_s     ! Total mass of PCM in LHS (kg)
        
    Half_air_gap = Air_gap/2.             ! Half thickness, which is what is taken in the model since we are considering only the heat transfer 
                                          ! to/from one slab and using symmetry (m)

    Dh = 4*Duct_width*Air_gap/(2*(Duct_width + Air_gap))    ! Hydraulic diameter of air channels, assuming the slabs go
                                                            ! all the way to the walls (m)
    
    ! CREATING THE NECESSARY VARIABLES, VECTORS AND ARRAYS FOR THE SCHEME
    M = nint(Slab_length/dx)              ! Number of divisions along slab length (x direction), -   
    N = nint(Slab_half_thickness/dy)      ! Number of divisions along half the slab thickness (y direction), - 
    dy = Slab_half_thickness/(N+0.5)

    allocate(x_position(0:M))               ! Defines the size of the vector x_position
    allocate(y_position(0:N))               ! Defines the size of the vector y_position
        
    do j=0,M                                ! Assign values of x_position according to slab length and M
        x_position(j) = j*Slab_length/M
    end do

    do i=0,N                                ! Assign values of y_position according to slab half thickness and N
        y_position(i) = dy/2+i*dy
    end do
    
    dx = x_position(1) - x_position(0)  ! Mesh size along slab lenth (actual) (m)
    dy = y_position(1) - y_position(0)  ! Mesh size along slab thickness (actual) (m)
    
    fe = (dy/2.)/(dy/2. + thickn_w)     ! Ratio of distances for harmonic mean interface conductivity between container 
                                        ! wall and PCM (Book by Patankar page. 45)
    
    allocate(Ta(0:M)); Ta = Ta0             ! Defines the size of the vector Ta and initializes to Ta0
    allocate(Ta_new(0:M)); Ta_new = Ta      ! Defines the size of the vector Ta_new and initializes to Ta0
    
    allocate(Tw(0:M)); Tw = Tw0             ! Defines the size of the vector Tw and initializes to Tw0
    allocate(Tw_new(0:M)); Tw_new = Tw      ! Defines the size of the vector Tw_new and initializes to Tw0
    
    allocate(T(0:N,0:M)); T = T0            ! Defines the size of the array T and initializes to T0
    allocate(T_new(0:N,0:M)); T_new = T     ! Defines the size of the array T_new and initializes to 0
    
    
    allocate(c(0:N,0:M))                    ! Matrix of specific heat values in the grid (J/kg-K)
    allocate(k(0:N,0:M))                    ! Matrix of thermal conductivity values in the grid (W/m-K)
    
    where(T >= Tm)                          ! Initializing values of matrices c and k 
        c = cl; k = kl
    elsewhere 
        c = cs;  k = ks
    endwhere
      
    allocate(H(0:N,0:M))                ! Defines the size of the array H 
    allocate(H_new(0:N,0:M))            ! Defines the size of the array H_new
    
    where(T <= Ts)                      ! Initializing values of matrix H 
        H = c*(T - Ts)
    elsewhere(T >= Tl)
        H = L + c*(T - Tl)
    elsewhere
        H = L*(T - Ts/(T_range))
    endwhere
    
    H_new = H                           ! Setting the matrix H_new equal to H
    
    allocate(Lambda(0:N,0:M))           ! Defines size of matrix to track liquid fraction
    
    where(H <= Hs)                      ! Defining the values of matrix Lamdba according to H
        Lambda = 0.0
    elsewhere(H >= Hl)
        Lambda = 1.0
    elsewhere
        Lambda = H/Hl
    endwhere   
    
    Lambda_avg = sum(Lambda)/size(Lambda)       ! Calculates the mean liquid fraction of the entire PCM region
    
    Store_dt = Store_dt_min*60                  ! Data storage interval (s)
    Store_dt_details = Store_dt_min_details*60  ! Data storage interval for detailed slab data (s)
    
    Total_sim_t = size(t_min)*Data_time_step    ! The total simulation time (s)
 
    allocate(delta_H(0:N,0:M)); delta_H = 1     ! Matrix of iteration changes of H in each node, inititalized to value above the tolerance desired
    allocate(delta_Ta(0:M)); delta_Ta = 1       ! Matrix of iteration changes of Ta in each node, initialized to value above the tolerance desired
    allocate(delta_Tw(0:M)); delta_Tw = 1       ! Matrix of iteration changes of Tw in each node, initialized to value above the tolerance desired
    
    iterations_tot = 0
    
    Curr_sim_time = 0.0                         ! The current simulation time during the calculations (s)
    
    Next_store_time = 0.          !Initializing the variable keeping track of the next simulation time to store values
    Next_store_time_details = 0.  !Initializing the variable keeping track of the next simulation time to store detail values
    
    ! HERE WE WRITE TO THE FILE THAT IS GOING TO BE USED BY THE DOWNSTREAM COMPONENT
    open(unit=output_unit_number, file=trim(output_file), status='replace', action='write')
    
    ! Write setup info to output files  
    write(output_unit_number, '(8(A,F10.2))') "rho_s:", rho_s, " rho_l:", rho_l, " ks:", ks, " kl:", kl, " cs:", cs, &
        " cl:", cl, " L:", L, " Tm:", Tm
    write(output_unit_number, '(A,I8,2(A,F8.2),3(A,F8.3),2(A,F8.1))') "N_slabs:", N_slabs, " Slab_Length:", Slab_length, &
        " Slab_width:", Slab_width, " Slab_thickness:", Slab_thickness, " PCM_mass:", mass_PCM, " Air_gap:", Air_gap, " T0:", T0, " Ta0:", Ta0
    write(output_unit_number, '(A,F8.1,A,F8.4)') "T_range:", T_range, " Wall_thick:", thickn_w
    write(output_unit_number, '(A,F8.1,2(A,F8.5))') "dt:", dt, " dx:", dx, " dy:", &
        dy
    write(output_unit_number, '(A12,A12,F9.3,F8.1,F8.1,F8.1,F8.1,F8.0,F8.4,F13.3,F11.3,F11.8)') Date(0), Hour_min(0), Hour_dec(0), Curr_sim_time/60., &
        Ins(0), Tamb(0), wind(0), Patm(0), X(0), Ta(M), Lambda_avg, mdot_a(0)

    write(output_unit_number, '(A,A12,A9,A8,A8,A8,A8,A8,A8,A13,A11,A11)') "Date", "Hour_min", "Hour_dec", "t_min", "Ins", "Tamb", "Wind", "Patm", "X", &
         "Tpcm_out_avg", "Lambda_avg", "mdot_a"
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! HERE WE WRITE TO THE FILE THAT IS GOING TO CONTAIN DETAILED TEMPERATURE IN SLAB
    open(unit=output_DETAILS_T_unit_number, file=trim(output_DETAILS_T_file), status='replace', action='write')
 
    write(output_DETAILS_T_unit_number, '(8(A,F10.2))') "rho_s:", rho_s, " rho_l:", rho_l, " ks:", ks, " kl:", kl, " cs:", cs, &
        " cl:", cl, " L:", L, " Tm:", Tm
    write(output_DETAILS_T_unit_number, '(A,I8,2(A,F8.2),3(A,F8.3),2(A,F8.1))') "N_slabs:", N_slabs, " Slab_Length:", Slab_length, &
        " Slab_width:", Slab_width, " Slab_thickness:", Slab_thickness, " PCM_mass:", mass_PCM, " Air_gap:", Air_gap, " T0:", T0, " Ta0:", Ta0
    write(output_DETAILS_T_unit_number, '(A,F8.1,A,F8.4,A,F8.1,2(A,F8.5))') "T_range:", T_range, " Wall_thick:", thickn_w, " dt:", dt, " dx:", dx, " dy:", &
        dy
    write(output_DETAILS_T_unit_number, '(A,F10.1)') "time,s:", Curr_sim_time
    write(output_DETAILS_T_unit_number, '(A)' ,advance = 'no'), "y,m/x,m"
    
    do j = 0,size(x_position)-1
        write(output_DETAILS_T_unit_number, '(F8.3)', advance = 'no') x_position(j)
    end do
    write(output_DETAILS_T_unit_number,*)
    
    write(output_DETAILS_T_unit_number, '(A)', advance = 'no'), "Ta"
    do j = 0,size(Ta_new)-1
        write(output_DETAILS_T_unit_number, '(F8.3)', advance = 'no') Ta(j)
    end do
    write(output_DETAILS_T_unit_number,*)
    
    write(output_DETAILS_T_unit_number, '(A)', advance = 'no'), "Tw"
    do j = 0,size(Tw_new)-1
        write(output_DETAILS_T_unit_number, '(F8.3)', advance = 'no') Tw(j)
    end do
    write(output_DETAILS_T_unit_number,*)
    
    do i = 0,N
        write(output_DETAILS_T_unit_number, '(F7.5)', advance = 'no') y_position(i)
        do j = 0,M
            write(output_DETAILS_T_unit_number, '(F8.3)', advance = 'no') T(i,j)
        end do
        write(output_DETAILS_T_unit_number,*)
    end do
    write(output_DETAILS_T_unit_number,*) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! HERE WE WRITE TO THE FILE THAT IS GOING TO CONTAIN DETAILED ENTHALPY IN SLAB
    open(unit=output_DETAILS_H_unit_number, file=trim(output_DETAILS_H_file), status='replace', action='write')
 
    write(output_DETAILS_H_unit_number, '(8(A,F10.2))') "rho_s:", rho_s, " rho_l:", rho_l, " ks:", ks, " kl:", kl, " cs:", cs, &
        " cl:", cl, " L:", L, " Tm:", Tm
    write(output_DETAILS_H_unit_number, '(A,I8,2(A,F8.2),3(A,F8.3),2(A,F8.1))') "N_slabs:", N_slabs, " Slab_Length:", Slab_length, &
        " Slab_width:", Slab_width, " Slab_thickness:", Slab_thickness, " PCM_mass:", mass_PCM, " Air_gap:", Air_gap, " T0:", T0, " Ta0:", Ta0
    write(output_DETAILS_H_unit_number, '(A,F8.1,A,F8.4,A,F8.1,2(A,F8.5))') "T_range:", T_range, " Wall_thick:", thickn_w, " dt:", dt, " dx:", dx, " dy:", &
        dy
    write(output_DETAILS_H_unit_number, '(A,F10.1)') "time,s:", Curr_sim_time
    write(output_DETAILS_H_unit_number, '(A)' ,advance = 'no'), "y,m/x,m"
    
    do j = 0,size(x_position)-1
        write(output_DETAILS_H_unit_number, '(F10.3)', advance = 'no') x_position(j)
    end do
    write(output_DETAILS_H_unit_number,*)
    
    write(output_DETAILS_H_unit_number, '(A)', advance = 'no'), "Ta"
    do j = 0,size(Ta_new)-1
        write(output_DETAILS_H_unit_number, '(F10.3)', advance = 'no') Ta(j)
    end do
    write(output_DETAILS_H_unit_number,*)
    
    write(output_DETAILS_H_unit_number, '(A)', advance = 'no'), "Tw"
    do j = 0,size(Tw_new)-1
        write(output_DETAILS_H_unit_number, '(F10.3)', advance = 'no') Tw(j)
    end do
    write(output_DETAILS_H_unit_number,*)
    
    do i = 0,N
        write(output_DETAILS_H_unit_number, '(F7.5)', advance = 'no') y_position(i)
        do j = 0,M
            write(output_DETAILS_H_unit_number, '(F10.1)', advance = 'no') H(i,j)
        end do
        write(13,*)
    end do
    write(13,*)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    A = 1/(rho_a*c_a*Half_air_gap)        ! Factor used in equation for the air to simplify writing
       
    Next_store_time = Next_store_time + Store_dt
    Next_store_time_details = Next_store_time_details + Store_dt_details
    
    Tout_avg_stor_dt = 0.0
    Count_stor_dt = 1           ! This counts how many timesteps have passed before storing results in the first output file, so that 
                                ! the avg of outlet Ta during the store_dtcan be calculated.
                                ! Starts at 1 because the least is to write data every time step.
    
    ! TIME LOOP STARTS HERE
    do while (Curr_sim_time < Total_sim_t)
        iterations = 0
        Tair_in = Ta_in(int(Curr_sim_time/Data_time_step))      ! Get current inlet air temperature
        mdot_air = mdot_a(int(Curr_sim_time/Data_time_step))    ! Get current air mass flow rate, kg/s
        V_a = mdot_air/rho_a                                    ! Airflow rate, m3/s
        va =  V_a/(N_channels*Air_gap*Duct_width)               ! Mean air velocity through gaps, m/s

        Re = rho_a*va*Dh/mu_a                                   ! Reynolds number in air channels 
        ff = (0.79*log(Re) - 1.64)**(-2)                        ! Friction factor to be used in Nu calculation with Gnielinski eq.
       
        ! Calculation of heat transfer coefficient between air and slab surface
        if (Re <= 2300.0) then
            Xaster = (Slab_length/Dh)/(Re*Pr)                       ! Factor for Nu in laminar flow 
            Nu = 7.55 + (0.024*Xaster**(-1.14))/(1 + 0.0358*(Pr**0.17)*Xaster**(-0.64)) ! Nu in parallel plates simultaneously developing flow ROHSENOW PAG 5.63
        else
            Nu = (ff/8.)*(Re-1000.)*Pr/(1 + 12.7*(ff/8.)**0.5*(Pr**(2./3.)-1.))     ! Nu by Gnielinski (appears in practically all sources as a good eq., 
                                                                                ! used by Dolado and others)
        end if

        h_aw = Nu*k_a/Dh
                
        do while (maxval(delta_H) > TOL_delta_H .or. maxval(delta_Ta) > TOL_delta_Ta .or. maxval(delta_Tw) > TOL_delta_Tw)  ! This is the iteration control, 
                                                                                                                            ! the tolerance value can be varied
            ! LOOP OVER AIR NODES FROM i=0 TO M 
            Ta_iter = Tair_in
            Ta_new(0) = Ta_iter
            delta_Ta(0) = Ta_iter - Ta_new(0)
            
            do j = 1,M
                Ta_iter = (Ta(j) + (va*dt/dx)*Ta_new(j-1) + dt*A*h_aw*Tw_new(j))/(1 + dt*(va/dx + A*h_aw))
                delta_Ta(j) = abs(Ta_iter - Ta_new(j))      ! Recording the change in temperature from past iteration
                Ta_new(j) = Ta_iter                         ! Updating value of Ta_new(j)
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

                Tw_iter = relax*Tw_iter + (1-relax)*Tw_new(j)             ! Applying underrelaxation 

                delta_Tw(j) = abs(Tw_iter - Tw_new(j))              ! Recording the change in temperature from past iteration
                
                Tw_new(j) = Tw_iter                         ! Updating value of Tw_new(j)
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
                
                kx_ant = 2*k(i,jm1)*k(i,j)/(k(i,jm1)+k(i,j))  !Patankar formulas for interface conductivity using harmonic mean
                kx_post = 2*k(i,jp1)*k(i,j)/(k(i,jp1)+k(i,j))
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
            
                
                H_new(i,j) = H_iter                 ! Updating value of H_new(i,j) and T_new(i,j)
                T_new(i,j) = T_iter
            
                if (H_new(i,j) <= Hs) then
                    k(i,j) = ks
                    c(i,j) = cs
                else if (Hs < H_new(i,j) .and. H_new(i,j) < Hl) then
                    k(i,j) = (kl - ks)/(Hl - Hs)*H_new(i,j) + ks
                    c(i,j) = (cl - cs)/(Hl - Hs)*H_new(i,j) + cs
                else
                    k(i,j) = kl
                    c(i,j) = cl
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
                    
                    kx_ant = 2*k(i,jm1)*k(i,j)/(k(i,jm1)+k(i,j))  ! Patankar formulas for interface conductivity using harmonic mean
                    kx_post = 2*k(i,jp1)*k(i,j)/(k(i,jp1)+k(i,j)) 
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

                    
                    H_new(i,j) = H_iter                      ! Updating value of H_new(i,j) and T_new(i,j)
                    T_new(i,j) = T_iter
                                    
                    if (H_new(i,j) <= Hs) then
                        k(i,j) = ks
                        c(i,j) = cs
                    else if (Hs < H_new(i,j) .and. H_new(i,j) < Hl) then
                        k(i,j) = (kl - ks)/(Hl - Hs)*H_new(i,j) + ks
                        c(i,j) = (cl - cs)/(Hl - Hs)*H_new(i,j) + cs
                    else
                        k(i,j) = kl
                        c(i,j) = cl
                    end if
                end do
            end do
            
            ! NOW THE SYMMETRY ROW i=N, GOING FROM COLUMN j=0 TO M 
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
                
                kx_ant = 2*k(i,jm1)*k(i,j)/(k(i,jm1)+k(i,j))   ! Patankar formulas for interface conductivity using harmonic mean
                kx_post = 2*k(i,jp1)*k(i,j)/(k(i,jp1)+k(i,j))  
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
                
                delta_H(i,j) = abs(H_iter - H_new(i,j))             ! Recording the change in enthalpy from past iteration
                
                H_new(i,j) = H_iter                             ! Updating value of H_new(i,j) and T_new(i,j)
                T_new(i,j) = T_iter
                
                if (H_new(i,j) <= Hs) then
                    k(i,j) = ks
                    c(i,j) = cs
                else if (Hs < H_new(i,j) .and. H_new(i,j) < Hl) then
                    k(i,j) = (kl - ks)/(Hl - Hs)*H_new(i,j) + ks
                    c(i,j) = (cl - cs)/(Hl - Hs)*H_new(i,j) + cs
                else
                    k(i,j) = kl
                    c(i,j) = cl
                end if
            end do
            
            iterations = iterations + 1
            iterations_tot = iterations_tot + 1
        
        end do
        
        print *, "max_delta_Ta", maxval(delta_Ta)
        print *, "max_delta_Tw", maxval(delta_Tw)
        print *, "max_delta_H", maxval(delta_H)
                
        print *, "t", Curr_sim_time
        print *, "iterations", iterations       
        
        delta_H = 1.0    ! Reset the values of the deltas in each node to a value above the tolerance desired
        delta_Ta = 1.0
        delta_Tw = 1.0        
        
        where(H_new <= Hs)                      ! Defining the values of matrix Lamdba according to H
            Lambda = 0.0
        elsewhere(H_new >= Hl)
            Lambda = 1.0
        elsewhere
            Lambda = H_new/Hl
        endwhere   
    
        Lambda_avg = sum(Lambda)/size(Lambda) ! Calculates the mean liquid fraction of the entire PCM region
        
        pcm_width = Slab_width + Lambda_avg*(Duct_width - Slab_width)   ! This PCM width is the actual "height" of the PCM inside
                                                    ! the slabs as it melts according to how much PCM has melted (liquid fraction).
                                                    ! It can go from the level when all PCM is solid, to the level when all PCM is 
                                                    ! liquid, according to the calculated average liquid fraction
        
        Ta_out_mixed = (Ta_new(M)*pcm_width + Tair_in*(Duct_width - pcm_width))/Duct_width    ! This calculates the mixed
                                                    ! temperature of the air at the outlet of the slab and the air which was
                                                    ! above the pcm level,which leaves de LHS as it arrived, with the inlet T
                                                    ! which, if if air comes from a solar air heater, is the SAH outlet temperature
        
        Tout_avg_stor_dt =  (Tout_avg_stor_dt*(Count_stor_dt - 1) + Ta_out_mixed)/Count_stor_dt    ! This is used to calculate the 
                                                    ! average outlet T during the store_dt. It calculates the mean of all Ta_out_mixed 
                                                    ! values between store_dts to give a value of mean Tout for the whole storage period. 
                                                    ! For example if i store data every minute and my ts is 10s, the code calculates 
                                                    ! a Tout every 10s, so here we calculate the mean of 6 Tout values which represent 
                                                    ! the mean Tout of that minute. 

        Curr_sim_time = Curr_sim_time + dt          ! Advance time step

        H = H_new                                   ! Replacing the new calculated enthalpy and temperatures into the arrays of current values to move into next timestep
        T = T_new
        Ta = Ta_new
        Tw = Tw_new

 
        ! STORING DATA IN OUTPUT FILES
        
        indx = int(Curr_sim_time/Data_time_step)
        
        if (mod(Curr_sim_time, Data_time_step) == 0) then
            indx = indx - 1
        end if    
            
        if (Next_store_time - Curr_sim_time < dt) then
            write(output_unit_number, '(A12,A12,F9.3,F8.1,F8.1,F8.1,F8.1,F8.0,F8.4,F13.3,F11.3,F11.8)') Date(indx), Hour_min(indx), Hour_dec(indx), &
                Curr_sim_time/60., Ins(indx), Tamb(indx), wind(indx), Patm(indx), X(indx), Tout_avg_stor_dt, Lambda_avg, mdot_a(indx)
            print *, "saving..."
            
            Tout_avg_stor_dt = 0
            Count_stor_dt = 0
            Next_store_time = Next_store_time + Store_dt
        end if
        
        if (Next_store_time_details - Curr_sim_time < dt) then           

            write(output_DETAILS_T_unit_number, '(A,F10.1)') "time,s:", Curr_sim_time
            write(output_DETAILS_T_unit_number, '(A)' ,advance = 'no'), "y,m/x,m"
    
            do j = 0,size(x_position)-1
                write(output_DETAILS_T_unit_number, '(F8.3)', advance = 'no') x_position(j)
            end do
            write(output_DETAILS_T_unit_number,*)
    
            write(output_DETAILS_T_unit_number, '(A)', advance = 'no'), "Ta"
            do j = 0,size(Ta_new)-1
                write(output_DETAILS_T_unit_number, '(F8.3)', advance = 'no') Ta_new(j)
            end do
            write(output_DETAILS_T_unit_number,*)
    
            write(output_DETAILS_T_unit_number, '(A)', advance = 'no'), "Tw"
            do j = 0,size(Tw_new)-1
                write(output_DETAILS_T_unit_number, '(F8.3)', advance = 'no') Tw_new(j)
            end do
            write(output_DETAILS_T_unit_number,*)
    
            do i = 0,N
                write(output_DETAILS_T_unit_number, '(F7.5)', advance = 'no') y_position(i)
                do j = 0,M
                    write(output_DETAILS_T_unit_number, '(F8.3)', advance = 'no') T_new(i,j)
                end do
                write(output_DETAILS_T_unit_number,*)
            end do
            write(output_DETAILS_T_unit_number,*) 
            
            write(output_DETAILS_H_unit_number, '(A,F10.1)') "time,s:", Curr_sim_time
            write(output_DETAILS_H_unit_number, '(A)' ,advance = 'no'), "y,m/x,m"
    
            do j = 0,size(x_position)-1
                write(output_DETAILS_H_unit_number, '(F10.3)', advance = 'no') x_position(j)
            end do
            write(output_DETAILS_H_unit_number,*)
    
            write(output_DETAILS_H_unit_number, '(A)', advance = 'no'), "Ta"
            do j = 0,size(Ta_new)-1
                write(output_DETAILS_H_unit_number, '(F10.3)', advance = 'no') Ta_new(j)
            end do
            write(output_DETAILS_H_unit_number,*)
    
            write(output_DETAILS_H_unit_number, '(A)', advance = 'no'), "Tw"
            do j = 0,size(Tw_new)-1
                write(output_DETAILS_H_unit_number, '(F10.3)', advance = 'no') Tw_new(j)
            end do
            write(output_DETAILS_H_unit_number,*)
    
            do i = 0,N
                write(output_DETAILS_H_unit_number, '(F7.5)', advance = 'no') y_position(i)
                do j = 0,M
                    write(output_DETAILS_H_unit_number, '(F10.1)', advance = 'no') H_new(i,j)
                end do
                write(output_DETAILS_H_unit_number,*)
            end do
            write(output_DETAILS_H_unit_number,*)
            
            
            Next_store_time_details = Next_store_time_details + Store_dt_details  
        end if        
        
        Count_stor_dt = Count_stor_dt + 1

    end do
deallocate(Date,Hour_min,Hour_dec,t_min,Ins,Tamb,wind,Patm,X,Ta_in,x_position,y_position,Ta,Ta_new)
deallocate(Tw,Tw_new,T,T_new,H,H_new,Lambda,c,k)
deallocate(delta_H,delta_Ta,delta_Tw)
    
end program LHS_2D