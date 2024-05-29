module file_data_module
    implicit none
    integer :: num_rows
    character(len=10), allocatable :: Date(:), Hour_min(:)
    double precision, allocatable :: Hour_dec(:), t_min(:), Ins(:), Tamb(:), wind(:), Patm(:), Xa_in(:), Ta_in(:), Lambda(:), mdot_a(:)
end module file_data_module

! HERE ARE TWO SUBROUTINES TO READ THE INPUT DATA FILE AND ASSIGN EACH DATA VARIABLE TO A LIST WHOSE LENGTH IS THE NUMBER OF DATA ROWS.
! AS IT IS, THE FIRST SUBROUTINE IS MADE TO READ THE OUPUT FILE FROM THE SOLAR AIR HEATER (SAH) PROGRAM WHICH HAS 6 ROWS WITH SETUP INFO, WHICH ARE SKIPPED:
! THE USER CAN MODIFY THE SUBROUTINE OR THE INPUT FILE TO MATCH IF THE DATA FILE FORMAT IS DIFFERENT
subroutine read_file_SAH_output(input_file)
    use file_data_module
    implicit none
    character(len=*), intent(in) :: input_file
    integer :: input_unit_number_SAH, i, ierr
    character(len=150) :: line

    input_unit_number_SAH = 10
    
    !Open the file
    open(unit=input_unit_number_SAH, file=input_file, status='old')

    !Skip the first 6 rows
    do i = 0, 5
        read(input_unit_number_SAH, *)
    end do
    
    !Initialize variables
    num_rows = 0
    ierr = 0

    !First pass: Count the number of rows
    do
        read(input_unit_number_SAH, '(A)', iostat=ierr)
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
    allocate(Xa_in(0:num_rows-1))
    allocate(Ta_in(0:num_rows-1))
    allocate(mdot_a(0:num_rows-1))
    
    !Rewind the file to the beginning
    rewind(input_unit_number_SAH)

    !Skip the first 6 rows again
    do i = 0, 5
        read(input_unit_number_SAH, *)
    end do
    
    !Second pass: Read data and store in arrays
    do i = 0, num_rows-1
        read(input_unit_number_SAH, '(A)') line
        !print *, line
        !Parse the line and store values in arrays
        read(line, *) Date(i), Hour_min(i), Hour_dec(i), t_min(i), Ins(i), Tamb(i), wind(i), Patm(i), Xa_in(i), Ta_in(i), mdot_a(i)
    end do

    !Close the file
    close(input_unit_number_SAH)
    
 end subroutine read_file_SAH_output
    
! THIS SECOND SUBROUTINE TO READ INPUT DATA IS MEANT TO USE THE OUTPUT DATA FROM THE LATENT HEAT STORAGE (LHS) AS INPUT. 
! IT DIFFERS FROM THE FIRST IN THAT THE LHS OUTPUT FILE HAS AN EXTRA COLUMN FOR THE LIQUID FRACTION WHICH IS NOT NEEDED HERE
! BUT IT HAS TO BE ACCOUNTED FOR    
subroutine read_file_LHS_output(input_file)
    use file_data_module
    implicit none
    character(len=*), intent(in) :: input_file
    integer :: input_unit_number_LHS, i, ierr
    character(len=150) :: line
    
    input_unit_number_LHS = 11
    
    ! Open the file
    open(unit=input_unit_number_LHS, file=input_file, status='old')

    ! Skip the first 6 rows
    do i = 0, 5
        read(input_unit_number_LHS, *)
    end do
    
    ! Initialize variables
    num_rows = 0
    ierr = 0

    ! First pass: Count the number of rows
    do
        read(input_unit_number_LHS, '(A)', iostat=ierr)
        if (ierr /= 0) exit  ! Exit loop if end-of-file or error
        num_rows = num_rows + 1
    end do

    ! Allocate memory for the arrays based on the number of rows
    allocate(Date(0:num_rows-1))
    allocate(Hour_min(0:num_rows-1))
    allocate(Hour_dec(0:num_rows-1))
    allocate(t_min(0:num_rows-1))
    allocate(Ins(0:num_rows-1))
    allocate(Tamb(0:num_rows-1))
    allocate(wind(0:num_rows-1))
    allocate(Patm(0:num_rows-1))
    allocate(Xa_in(0:num_rows-1))
    allocate(Ta_in(0:num_rows-1))
    allocate(Lambda(0:num_rows-1))
    allocate(mdot_a(0:num_rows-1))

    ! Rewind the file to the beginning
    rewind(input_unit_number_LHS)

    ! Skip the first 6 rows again
    do i = 0, 5
        read(input_unit_number_LHS, *)
    end do

    ! Second pass: Read data and store in arrays
    do i = 0, num_rows-1
        read(input_unit_number_LHS, '(A)') line
        !print *, line
        ! Parse the line and store values in arrays
        read(line, *) Date(i), Hour_min(i), Hour_dec(i), t_min(i), Ins(i), Tamb(i), wind(i), Patm(i), Xa_in(i), Ta_in(i), Lambda(i), mdot_a(i)
    end do

    ! Close the file
    close(input_unit_number_LHS)    
 end subroutine read_file_LHS_output

! FUNCTION TO CALCULATE THE SATURATION VAPOR PRESSURE FROM A GIVEN TEMPERATURE 
double precision function saturation_vapor_pressure(Tad)
    implicit none
    double precision :: Tad, Pvsat, Tabs, log_Pvsat
    
    Tabs = Tad + 273.15
    log_Pvsat = -5.8002206e3/Tabs + 1.3914493e0 - 4.8640239e-2*Tabs + 4.1764768e-5*Tabs**2 - 1.4452093e-8*Tabs**3 + 6.5459673e0*log(Tabs)
    saturation_vapor_pressure = exp(log_Pvsat)
    return
end function saturation_vapor_pressure    

    
! HERE STARTS THE MAIN PROGRAM TO CALCULATE THE OPERATION OF THE DEEP BED DRYER    
program Deep_Bed_Drying
    use file_data_module
    implicit none
    
    double precision :: relax, hfg, Xp0_wetting
    double precision :: C1, C2, C3 ! Oswin isortherm parameters
    double precision :: Mass, Xp0, porosity, ap, Dp, Side_length
    double precision :: rho_a, c_a, c_v, c_w, mu_a, dz, dt, p_atm, rho_da
    double precision :: TOL_delta_Ta, TOL_delta_Xa, TOL_delta_Tp, TOL_delta_Xp
    double precision :: rho_bulk, rho_dp_bulk, Bulk_volume, c_dp, K, N
    double precision :: Ta0, Tp0, Mp0, Ta_iter, Xa_iter, Tp_iter, Xp_iter, Xp_eq, Xp_current, Xp_init, MR, teq
    double precision :: pvsat_Ta, pv_Ta, pv_Ta_in, pvsat_Ta_in, rh, rh0, pvsat0, pv0, Xa0, A, B, C, saturation_vapor_pressure   
    double precision :: Height, Area_bin, Vdot_a, v_a, Re, h, Tair_in, mdot_air, Xair_in
    double precision :: Curr_sim_time, Data_time_step, Store_dt, Store_dt_min, Total_sim_t, Next_store_time
    integer :: i, j, ind, index, iterations, iterations_tot, numb_stor_div, M, output_unit_number
    double precision, allocatable :: rh_in(:), z_position(:)
    double precision, allocatable :: Ta(:), Ta_new(:), Xa(:), Xa_new(:), Tp(:), Tp_new(:), Xp(:), Xp_new(:)
    integer, allocatable :: write_pos_index(:)
    double precision, allocatable :: delta_Ta(:), delta_Xa(:), delta_Tp(:), delta_Xp(:)
    
    character(len=100) :: output_file = 'Output_Dryer_from_LHS_Stutt_Aug2022.txt'
    
    output_unit_number = 12
    
    Data_time_step = 60.    ! This is the time resolution of the input data in seconds
    relax = 0.9             ! Underrelaxation factor, can be lower than 1 for stability if required but more than 0
    
    !Call the subroutine to read data
    call read_file_LHS_output("LHS_output_Stutt_Aug2022_as_INPUT.txt")
    
    ! Make a rh_inlet vector based on the input data air conditions
    allocate(rh_in(0:size(Xa_in)-1))
    do i = 0,(size(Xa_in)-1)
        pv_Ta_in = Xa_in(i)*Patm(i)/(0.621945 + Xa_in(i))
        pvsat_Ta_in = saturation_vapor_pressure(Ta_in(i))
        rh_in(i) = pv_Ta_in/pvsat_Ta_in
    end do
    
    hfg = 2253430.              ! Latent heat of vaporization (assumed constant) (J/kg)
    Xp0_wetting = 0.0            

    C1 = 0.129                  ! Parameters of Oswin isotherm equation for pionier WHEAT (paper Ramaj, 2021)
    C2 = -6.46e-4               ! This has to be changed by the user to parameters of the equation for his/her own product (there are various isotherm model equations)
    C3 = 2.944                  ! Care has to be taken not only with which isotherm equation is being used 
                                !  but also the units of the variables (rh as % or decimal, moisture as % or decimal and as wet basis or dry basis)
    
    ! INITIAL CONDITIONS
    Mass = 1000.                 ! Mass of product at initial moisture content (kg)
    Xp0 = 0.28205                 ! Initial product moisture content, db, (kg/kg)
    Tp0 = Tamb(0)                   ! Initial product temperature ( C)
    Ta0 = Tamb(0)                  ! Initial interstitial air temperature ( C)

    !  PRODUCT PROPERTIES
    porosity = 0.4                                 ! Porosity of the product bulk, from ASAE D241.4, -
    rho_bulk = 774.4 - 703*(Xp0/(1+Xp0)) + 18510*(Xp0/(1+Xp0))**2 - 148960*(Xp0/(1+Xp0))**3 + 311600*(Xp0/(1+Xp0))**4   ! Product bulk density depending on initial moisture (Nelson equation in ASAE D241.4 ) (kg/m3)
                                                                                                                        ! This can be replaced by a single value
    
    ap = 1181.          ! Product specific volumetric area of product (Brooker, Bakker A. Table 2.6, 8.1 and D-5), m2/m3
    c_dp = 1240.        ! Product specific heat of product in dry basis, from ASAE D243.4, J/kg-K
    Dp = 0.004          ! Product equivalent particle diameter for Reynolds number for heat trans. coeff. (from Brooker Tables 2.1 and 8.1), m
    
    ! AIR PROPERTIES
    rho_a = 1.127               ! Air density (at about 40  C) (kg/m3)
    c_a = 1007.                 ! Specific heat of dry air (J/kg-K)
    c_v = 1883.                 ! Specific heat of water vapor, (J/kg-K)
    mu_a = 1.907e-5             ! Air viscosity (Pa s) 
    
    c_w = 4187.                 ! Specific heat of water (J/kg-K)
    
    !  DRYER DIMENSIONS
    Side_length = 1.5           ! Assumed that the dryer is of square shape (m)
    Area_bin = Side_length**2   ! Bin floor area (m2) (if rectangular bin)
    !Dryer_diameter = 1.2       ! Dryer diameter (m) (if round bin)
    !Area_bin = pi*(Dryer_diameter/2)**2 ! Bin floor area (m2) (if round bin)
    
    ! SPACE AND TIME GRID SIZE
    dt = 30.                          ! Timestep (s)
    dz = 0.008                       ! Mesh size along dryer height (provisional) (m)
    
    ! CHOOSE HOW OFTEN TO STORE DATA
    Store_dt_min = 1.0               ! Defines how often data will be saved to file (min)
    
    ! TOLERANCE FOR THE DIFFERENCE BETWEEN CONSECUTIVE ITERATIONS
    TOL_delta_Ta = 0.001        ! Tolerance of nodal air temperature ( C)
    TOL_delta_Xa = 0.00001      ! Tolerance of nodal air humidity ratio (kg/kg)
    TOL_delta_Tp = 0.001        ! Tolerance of nodal product temperature ( C)
    TOL_delta_Xp = 0.00001      ! Tolerance of nodal product moisture content (kg/kg)
    
    
    ! CALCULATIONS START HERE
    
    Mp0 = Xp0/(1+Xp0)           ! Initial product moisture content, wb (kg/kg)

    rh0 = (1./(((C1 + C2*Ta0)/Xp0)**C3 + 1))        ! Relative humidity of interstitial air in equilibrium with product at beginning, -
    pvsat0 = saturation_vapor_pressure(Ta0)         ! Saturation vapor pressure at the initial interstitial temperature (Pa)
    pv0 = rh0*pvsat0                                ! Actual initial vapor pressure of the interstitial air (Pa)
    Xa0 = 0.621945*pv0/(Patm(0) - pv0)              ! Initial humidity ratio of interstitial air (kg/kg)

    rho_dp_bulk = rho_bulk*(1-Mp0)               ! WHEAT bulk dry matter density (kg/m3)
                                               ! Since the volume of the grain bed is fixed in the model the dry matter bulk density is constant
    Bulk_volume = Mass/rho_bulk                 ! Volume of bulk (m3) 
                                    
    Height = Bulk_volume/Area_bin  ! Height of the grain bed (m) 
                                                                    
    ! CREATING THE NECESSARY VARIABLES, VECTORS AND ARRAYS FOR THE SCHEME
    M = int(Height/dz)          ! Number of divisions along dryer height
    numb_stor_div = 8             ! This is the number of divisions along the bed height in which i want to store data, first position being bottom and last being top of the bed
    
    do while (MOD(M,numb_stor_div) .NE. 0)  ! This is to find a number of divisions for the bed that is multiple of numb_stor_div
        M = M + 1
    end do
    
    dz = Height/M              ! Mesh size (actual), m
    allocate(z_position(0:M))   ! Defines the size of the vector z_position
    
    do i=0,M
        z_position(i) = i*Height/M ! Calculates node coordinates along the dryer height (m)
    end do
    
    allocate(write_pos_index(0:numb_stor_div))
    do i=0,numb_stor_div                            ! Here a list of indexes of the z_position vector is created where we want to store results
        write_pos_index(i) = int(i*1.0/(numb_stor_div)*M)
    end do
    
    allocate(Ta(0:M)); Ta = Ta0             ! Vector of interstitial air temperature along dryer height ( C)
    allocate(Ta_new(0:M)); Ta_new = Ta0     ! Vector of interstitial air temperature along dryer height at new time level ( C)
    allocate(Xa(0:M)); Xa = Xa0             ! Vector of interstitial air humidity along dryer height ( C)
    allocate(Xa_new(0:M)); Xa_new = Xa0     ! Vector of interstitial air humidity along dryer height at new time level ( C)
    allocate(Tp(0:M)); Tp = Tp0             ! Vector of product temperature along dryer height ( C)
    allocate(Tp_new(0:M)); Tp_new = Tp0     ! Vector of product temperature along dryer height in the next time level ( C)
    allocate(Xp(0:M)); Xp = Xp0             ! Vector of product moisture along dryer height ( C)
    allocate(Xp_new(0:M)); Xp_new = Xp0     ! Vector of product moisture along dryer height in the next time level ( C)
    
    
    
    allocate(delta_Ta(0:M)); delta_Ta = 1               ! Matrix of iteration changes of Ta in each node, inititalized to a value above the tolerance desired
    allocate(delta_Xa(0:M)); delta_Xa = 1               ! Matrix of iteration changes of Xa in each node, initialized to a value above the tolerance desired
    allocate(delta_Tp(0:M)); delta_Tp = 1               ! Matrix of iteration changes of Tp in each node, inititalized to a value above the tolerance desired
    allocate(delta_Xp(0:M)); delta_Xp = 1               ! Matrix of iteration changes of Xp in each node, initialized to a value above the tolerance desired
    
    Store_dt = Store_dt_min*60  ! Data storage interval in seconds (s)
    Total_sim_t = size(t_min)*Data_time_step      ! The total simulation time (s)
    
    iterations_tot = 0
    Curr_sim_time = 0
    Next_store_time = 0.
    
    ! WRITE SETUP INFO AND OUTPUT DATA TO FILE
    open(unit=output_unit_number, file=trim(output_file), status='replace', action='write')
    write(output_unit_number, '(6(A,F10.2))') "Mass:", Mass, " Volume:", Bulk_volume, " Xp0", Xp0, " Mp0:", Mp0, " Tp0:", Tp0, " porosity:", porosity
    write(output_unit_number, '(5(A,F8.2))') "rho_bulk:", rho_bulk, " rho_dp_bulk:", rho_dp_bulk, &
        " ap:", ap, " c_dp:", c_dp, " Dp:", Dp
    write(output_unit_number, '(5(A,F8.2))') "Side_length:", Side_length, " Area_bin:", Area_bin, " Height:", Height, " Ta0:", Ta0, " Xa0:", &
        Xa0
    write(output_unit_number, '(2(A,F8.2))') " dt:", dt, " dz:", dz
    
    write(output_unit_number, '(A,A12,A9,A8,A8,A8,A8,A8,A8,A8,A10)', advance ='no') "Date", "Hour_min", "Hour_dec", "t_min", "Ins", "Tamb", "Wind", "Patm", "Xa_in", &
         "Ta_in", "rh_in"
           
    do ind = 0,size(write_pos_index)-1
        write(output_unit_number, '(F8.4)', advance = 'no') z_position(write_pos_index(ind))
    end do
    write(output_unit_number, '(A)', advance = 'no') " Xp_avg"
    
    do ind = 0,size(write_pos_index)-1
        write(output_unit_number, '(F8.4)', advance = 'no') z_position(write_pos_index(ind))
    end do
    write(output_unit_number, '(A)', advance = 'no') " Tp_avg"
    
    do ind = 0,size(write_pos_index)-1
        write(output_unit_number, '(F8.4)', advance = 'no') z_position(write_pos_index(ind))
    end do    
    write(output_unit_number, '(A)', advance = 'no') " Xa_avg"
    
    do ind = 0,size(write_pos_index)-1
        write(output_unit_number, '(F8.4)', advance = 'no') z_position(write_pos_index(ind))
    end do
    write(output_unit_number, '(A)', advance = 'no') " Ta_avg"
    write(output_unit_number,*)
    
    write(output_unit_number, '(A,A12,F9.3,F8.1,F8.1,F8.1,F8.1,F8.0,F8.4,F8.3,F10.3)', advance ='no') Date(0), Hour_min(0), Hour_dec(0), Curr_sim_time/60., &
        Ins(0), Tamb(0), wind(0), Patm(0), Xa_in(0), Ta_in(0), rh_in(0)
    
    do ind = 0,size(write_pos_index)-1
        write(output_unit_number, '(F8.4)', advance = 'no') Xp(write_pos_index(ind))
    end do
    write(output_unit_number, '(F8.3)', advance = 'no') sum(Xp)/size(Xp)

    do ind = 0,size(write_pos_index)-1
        write(output_unit_number, '(F8.2)', advance = 'no') Tp(write_pos_index(ind))
    end do
    write(output_unit_number, '(F8.2)', advance = 'no') sum(Tp)/size(Tp)
    
    do ind = 0,size(write_pos_index)-1
        write(output_unit_number, '(F8.4)', advance = 'no') Xa(write_pos_index(ind))
    end do
    write(output_unit_number, '(F8.4)', advance = 'no') sum(Xa)/size(Xa)

    do ind = 0,size(write_pos_index)-1
        write(output_unit_number, '(F8.2)', advance = 'no') Ta(write_pos_index(ind))
    end do
    write(output_unit_number, '(F8.2)', advance = 'no') sum(Ta)/size(Ta)
    write(output_unit_number,*)
    
        
    Next_store_time = Next_store_time + Store_dt        ! Advance storage time

    
    ! TIME LOOP STARTS HERE
    do while (Curr_sim_time < Total_sim_t)
        iterations = 0
        Tair_in = Ta_in(int(Curr_sim_time/Data_time_step))      ! Get current inlet air temperature ( C)
        Xair_in = Xa_in(int(Curr_sim_time/Data_time_step))      ! Get current humidity ratio of inlet air (kg/kg)
        p_atm = Patm(int(Curr_sim_time/Data_time_step))         ! Get current atmospheric pressure (Pa)
        mdot_air = mdot_a(int(Curr_sim_time/Data_time_step))    ! Get current air mass flow rate (kg/s)
        
        Vdot_a = mdot_air/rho_a                                 ! Airflow rate, m3/s
        v_a = Vdot_a/Area_bin                                   ! Superficial velocity of drying air, m/s
        Re = v_a*rho_a*Dp/mu_a                                  ! Reynolds number
        h = 0.2755*(c_a/1000.)*(v_a*3600)*rho_a*Re**(-0.34)     ! Heat transfer coefficient according to Brooker, B-A book and paper de B-A 1967
        A = dt*v_a/(porosity*dz)                                ! Repeated factor in eq for Ta, Xa
        
        do while (maxval(delta_Ta) > TOL_delta_Ta .or. maxval(delta_Xa) > TOL_delta_Xa .or. maxval(delta_Tp) > TOL_delta_Tp .or. maxval(delta_Xp) > TOL_delta_Xp)
            ! LOOP OVER Ta NODES FROM i=0 TO M
            Ta_iter = Tair_in
            Ta_new(0) = Ta_iter
            delta_Ta(0) = Ta_iter - Ta_new(0)
            
            do j = 1,M
                rho_da = 1./(0.287042*(Ta_new(j) + 273.15)*(1 + 1.607858*Xa_new(j))/(0.001*p_atm)) ! Calculate dry air density for current node. First trying with Ta[j] and Xa[j] instead of new values.
                             
                B = porosity*rho_da*(c_a + c_v*Xa_new(j)) ! Repeated denominator of eq for Ta
                                            
                Ta_iter = (Ta(j) + A*Ta_new(j-1) + dt*h*ap*Tp_new(j)/B)/(1. + A + dt*h*ap/B)   
                                                                                                                                                                      
                Ta_iter = relax*Ta_iter + (1-relax)*Ta_new(j)       ! Applying underrelaxation 
                delta_Ta(j) = abs(Ta_iter - Ta_new(j))              ! Recording the change in temperature from past iteration
                Ta_new(j) = Ta_iter            
            end do
            
            ! LOOP OVER Xa NODES FROM i=0 TO M
            Xa_iter = Xair_in
            Xa_new(0) = Xa_iter
            delta_Xa(0) = Xa_iter - Xa_new(0)
        
            do j = 1,M  
                rho_da = 1/(0.287042*(Ta_new(j) + 273.15)*(1 + 1.607858*Xa_new(j))/(0.001*p_atm)) ! Calculate dry air density for current node
             
                
                Xa_iter = (Xa(j) + A*Xa_new(j-1) - rho_dp_bulk/(porosity*rho_da)*(Xp_new(j) - Xp(j)))/(1 + A)   ! Here also using Xp upwind values since air humidity will be
                                                                                                  ! more influenced by the product moisture content upwind of the air
                                                                                                                 ! node than by the downwind product node
                Xa_iter = relax*Xa_iter + (1-relax)*Xa_new(j)       ! Applying underrelaxation 
                delta_Xa(j) = abs(Xa_iter - Xa_new(j))              ! Recording the change in humidity ratio from past iteration
                Xa_new(j) = Xa_iter
            end do
            
            ! LOOP OVER Tp NODES FROM i=0 TO M 
            j = 0
            rho_da = 1/(0.287042*(Ta_new(j) + 273.15)*(1 + 1.607858*Xa_new(j))/(0.001*p_atm)) ! Calculate dry air density for current node
            
            C = (1 + (h*ap*dt*dz - rho_da*v_a*c_v*dt*(Xa_new(j+1) - Xa_new(j)))/(rho_dp_bulk*dz*(c_dp + c_w*Xp_new(j))))
          
            Tp_iter = (Tp(j) + (dt/(rho_dp_bulk*(c_dp + c_w*Xp_new(j))))*(h*ap*Ta_new(j) - (hfg + c_v*Ta_new(j))*(rho_da*v_a*(Xa_new(j+1) - Xa_new(j))/dz)))/C
               
            Tp_iter = relax*Tp_iter + (1-relax)*Tp_new(j)           ! Applying underrelaxation  
            delta_Tp(j) = abs(Tp_iter - Tp_new(j))                  ! Recording the change in temperature from past iteration
            Tp_new(j) = Tp_iter
            
            do j = 1,M-1
            
                rho_da = 1/(0.287042*(Ta_new(j) + 273.15)*(1 + 1.607858*Xa_new(j))/(0.001*p_atm)) ! Calculate dry air density for current node
            
                C = (1 + (2*h*ap*dt*dz - rho_da*v_a*c_v*dt*(Xa_new(j+1) - Xa_new(j-1)))/(2*rho_dp_bulk*dz*(c_dp + c_w*Xp_new(j))))
              
                Tp_iter = (Tp(j) + (dt/(rho_dp_bulk*(c_dp + c_w*Xp_new(j))))*(h*ap*Ta_new(j) - (hfg + c_v*Ta_new(j))*(rho_da*v_a*(Xa_new(j+1) - Xa_new(j-1))/(2*dz))))/C
                
                Tp_iter = relax*Tp_iter + (1-relax)*Tp_new(j)       ! Applying underrelaxation    
                delta_Tp(j) = abs(Tp_iter - Tp_new(j))              ! Recording the change in temperature from past iteration
                Tp_new(j) = Tp_iter
            end do
            
            j = M
            rho_da = 1/(0.287042*(Ta_new(j) + 273.15)*(1 + 1.607858*Xa_new(j))/(0.001*p_atm)) ! Calculate dry air density for current node
            
            C = (1 + (h*ap*dt*dz - rho_da*v_a*c_v*dt*(Xa_new(j) - Xa_new(j-1)))/(rho_dp_bulk*dz*(c_dp + c_w*Xp_new(j))))

            Tp_iter = (Tp(j) + (dt/(rho_dp_bulk*(c_dp + c_w*Xp_new(j))))*(h*ap*Ta_new(j) - (hfg + c_v*Ta_new(j))*(rho_da*v_a*(Xa_new(j) - Xa_new(j-1))/dz)))/C
        
            Tp_iter = relax*Tp_iter + (1-relax)*Tp_new(j)           ! Applying underrelaxation
        
            delta_Tp(j) = abs(Tp_iter - Tp_new(j))                  ! Recording the change in temperature from past iteration
        
            Tp_new(j) = Tp_iter  
            
            ! LOOP OVER PRODUCT MOISTURE CONTENT FROM i=0 TO M    
            ! HERE IT IS MADE USE OF A THIN LAYER DRYING EQUATION WHICH IS PRODUCT SPECIFIC
            ! AS WITH SORPTION ISOTHERMS, THERE ARE NUMEROUS MODELS. ONE USUALLY ADEQUATE FOR AGRIC. PRODUCTS IS THE PAGE MODEL; WHICH IS USED BELOW
            ! IF THE SIMPLRE LEWIS MODEL IS TO BE USED, THE PARAMETER N SHOULD BE SET TO 1.
            ! AS WITH THE ISOTHERM MODEL, CARE MUST BE TAKEN TO ENSURE THAT THE UNITS ARE TRANSFORMED TO THE UNITS OF THE PROGRAM IF REQUIRED.
            ! MOISTURE CONTENT SHOULD BE IN DRY BASIS, TIME IN SECONDS
            ! AS IT IS, THE PROGRAM USES THE SAME DRYING EQUATION FOR REWETTING. IF A REWETTING EQUATION IS APPROPRIATE THE USER CAN CAREFULLY MODIFY THE CODE
            
            do j = 0,M
            
                pvsat_Ta = saturation_vapor_pressure(Ta_new(j))  ! Saturation vapor pressure at the air temperature in the current node (Pa)
                                
                pv_Ta = Xa_new(j)*p_atm/(0.621945 + Xa_new(j))          ! Actual vapor pressure of air in current node (Pa)
                
                rh = pv_Ta/pvsat_Ta 
                              
                if (rh >= 0.99) then
                    rh = 0.99
                end if
                !Xp_eq = 0.01*(-log(1. - rh)/(AH*(Ta_new[j] + CH)))**(1./BH)
                Xp_eq = (C1 + C2*Ta_new(j))*(rh/(1.-rh))**(1./C3)
                !if Xp_eq != Xp_eq:
                 !   sys.exit()
                
                K = 2.8e-3*exp(0.059*Ta_new(j))*(100*rh)**(-0.139)*v_a**0.025   !Thin layer Page drying equation for wheat from Ramaj, 2021
                N = 0.784
            
            
                Xp_current = Xp_new(j)
                Xp_init = Xp0
                !if Xp_current > Xp0:
                 !   Xp_current = Xp0
            
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
                !teq = Curr_sim_time
                
                Xp_iter = (Xp(j) + (dt/60.)*N*K*((teq+dt)/60.)**(N-1)*Xp_eq)/(1 + (dt/60.)*N*K*((teq+dt)/60.)**(N-1))
                
                Xp_iter = relax*Xp_iter + (1-relax)*Xp_new(j)           ! Applying underrelaxation  
            
                delta_Xp(j) = abs(Xp_iter - Xp_new(j))                  ! Recording the change in moisture content from past iteration
            
            
                Xp_new(j) = Xp_iter
            
            end do
                
            iterations = iterations + 1
            iterations_tot = iterations_tot + 1
       
        end do 
            
        print *, "max_delta_Ta", maxval(delta_Ta)
        print *, "max_delta_Xa", maxval(delta_Xa)
        print *, "max_delta_Tp", maxval(delta_Tp)
        print *, "max_delta_Xp", maxval(delta_Xp)
        print *, "t", Curr_sim_time
        print *, "iterations", iterations
    
        delta_Ta = 1.0              ! Reset the values of the deltas in each node to a value above the tolerance desired
        delta_Xa = 1.0
        delta_Tp = 1.0
        delta_Xp = 1.0
    
        Curr_sim_time = Curr_sim_time + dt      ! Advance time step
    
        Ta = Ta_new                 ! Replacing the new calculated values into the vectors of current temperature to move into next timestep
        Xa = Xa_new
        Tp = Tp_new
        Xp = Xp_new
        
        
        ! STORING DATA IN OUTPUT FILE
        index = int(Curr_sim_time/Data_time_step)
        if (MOD(Curr_sim_time,Data_time_step) == 0) then
            index = index - 1
        end if
        
        if ((Next_store_time - Curr_sim_time) < dt) then
            write(output_unit_number, '(A,A12,F9.3,F8.1,F8.1,F8.1,F8.1,F8.0,F8.4,F8.3,F10.3)', advance ='no') Date(index), Hour_min(index), Hour_dec(index), Curr_sim_time/60., Ins(index), Tamb(index), wind(index), Patm(index), Xa_in(index), Ta_in(index), rh_in(index)
    
            do ind = 0,size(write_pos_index)-1
                write(output_unit_number, '(F8.4)', advance = 'no') Xp(write_pos_index(ind))
            end do
            write(output_unit_number, '(F8.3)', advance = 'no') sum(Xp)/size(Xp)

            do ind = 0,size(write_pos_index)-1
                write(output_unit_number, '(F8.2)', advance = 'no') Tp(write_pos_index(ind))
            end do
            write(output_unit_number, '(F8.2)', advance = 'no') sum(Tp)/size(Tp)
    
            do ind = 0,size(write_pos_index)-1
                write(output_unit_number, '(F8.4)', advance = 'no') Xa(write_pos_index(ind))
            end do
            write(output_unit_number, '(F8.4)', advance = 'no') sum(Xa)/size(Xa)

            do ind = 0,size(write_pos_index)-1
                write(output_unit_number, '(F8.2)', advance = 'no') Ta(write_pos_index(ind))
            end do
            write(output_unit_number, '(F8.2)', advance = 'no') sum(Ta)/size(Ta)
            write(output_unit_number,*)
        
            Next_store_time = Next_store_time + Store_dt
        end if
    end do
        
    deallocate(rh_in,z_position)
    deallocate(Ta,Ta_new,Xa,Xa_new,Tp,Tp_new,Xp,Xp_new)
    deallocate(write_pos_index)
    deallocate(delta_Ta,delta_Xa,delta_Tp,delta_Xp)
    deallocate(Date,Hour_min,Hour_dec,t_min,Ins,Tamb,wind,Patm,Xa_in,Ta_in)
    
        
end program Deep_Bed_Drying