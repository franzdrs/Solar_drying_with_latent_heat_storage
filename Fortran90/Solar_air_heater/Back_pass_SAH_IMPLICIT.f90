module file_data_module
    implicit none
    integer :: num_rows
    character(len=10), allocatable :: Date(:), Hour_min(:)
    double precision, allocatable :: Hour_dec(:), time_row(:), Ins(:), Tamb(:), rh_amb(:), wind(:), Patm(:), Xamb(:), Tdewp(:), mdot_a(:)
end module file_data_module

! FUNCTION TO CALCULATE THE SATURATION VAPOR PRESSURE AT A GIVEN TEMPERATURE
double precision function saturation_vapor_pressure(Tad)
    implicit none
    double precision :: Tad, Pvsat, Tabs, log_Pvsat
    
    Tabs = Tad + 273.15
    log_Pvsat = -5.8002206e3/Tabs + 1.3914493e0 - 4.8640239e-2*Tabs + 4.1764768e-5*Tabs**2 - 1.4452093e-8*Tabs**3 + 6.5459673e0*log(Tabs)
    saturation_vapor_pressure = exp(log_Pvsat)
    return
end function saturation_vapor_pressure       

! SUBROUTINE TO READ THE WEATHER DATA FILE AND ASSIGN EACH DATA VARIABLE TO A LIST WHOSE LENGTH IS THE NUMBER OF DATA ROWS.
! THE USER CAN PASS THE START ROW AND END ROW TO BBE READ FROM THE DATA FILE
! THE USER CAN MODIFY THE SUBROUTINE OR THE INPUT FILE TO MATCH IF THE INPUT DATA FILE FORMAT IS DIFFERENT    
subroutine read_file(input_file, start_row, end_row)
    use file_data_module
    implicit none
    character(len=*), intent(in) :: input_file
    integer :: input_unit_number, i, ierr, start_row, end_row, number_rows_to_read
    double precision :: Pvsat, Pv, saturation_vapor_pressure
    character(len=150) :: line

    input_unit_number = 10
    number_rows_to_read = end_row - start_row + 1
    ! Open the file
    open(unit=input_unit_number, file=input_file, status='old')

    ! Skip rows up to the start_row
    do i = 1,start_row-1    
        read(input_unit_number, *)
    end do
    
    ! Allocate memory for the arrays based on the number of rows
    allocate(Date(0:number_rows_to_read - 1))
    allocate(Hour_min(0:number_rows_to_read - 1))
    allocate(Hour_dec(0:number_rows_to_read - 1))
    allocate(time_row(0:number_rows_to_read - 1))
    allocate(Ins(0:number_rows_to_read - 1))
    allocate(Tamb(0:number_rows_to_read - 1))
    allocate(rh_amb(0:number_rows_to_read - 1))
    allocate(wind(0:number_rows_to_read - 1))
    allocate(Patm(0:number_rows_to_read - 1))
    allocate(Xamb(0:number_rows_to_read - 1))
    allocate(Tdewp(0:number_rows_to_read - 1))
    allocate(mdot_a(0:number_rows_to_read - 1))
        
    ! Read data and store in arrays
    do i = 0, number_rows_to_read - 1
        read(input_unit_number, '(A)') line
        ! Parse the line and store values in arrays
        read(line, *) Date(i), Hour_min(i), Hour_dec(i), time_row(i), Ins(i), Tamb(i), rh_amb(i), wind(i), Patm(i), mdot_a(i)
    end do
     
    do i = 0, number_rows_to_read - 1
        Pvsat = saturation_vapor_pressure(Tamb(i))
        Pv = (rh_amb(i)/100.)*Pvsat
        Xamb(i) = 0.621945*Pv/(Patm(i) - Pv)
        
        Tdewp(i) = (6.54 + 14.526*log(Pv/1000) + 0.7389*log(Pv/1000)**2 + 0.09486*log(Pv/1000)**3 + 0.4569*(Pv/1000)**0.1984) ! from ASHRAE (2017) Fundamentals pyschrometrics, eq 37
    end do
    ! Close the file
    close(input_unit_number)    
end subroutine read_file

    
! MAIN PROGRAM    
program SAH
    use file_data_module
    implicit none
    
    double precision :: sigma, pi, relax
    double precision :: h_dec, I, Ta, v_wind, p_atm, Xa, Tdp, Ts
    double precision :: L, W, H, gap_gp, tilt, Dh
    double precision :: thick_g, thick_p, thick_b, k_g, k_p, k_b, c_g, c_p, c_b, rho_g, rho_p, rho_b, Ub
    double precision :: thick_i, k_i, alpha_g, alpha_p, alpha_b, tau_g, e_p_up, e_p_down, e_g, e_b
    double precision :: rho_f, c_f, mu_f, k_f, Pr, th_diff_f
    double precision :: mdot_f, Vdot
    double precision :: Ac, v, spec_Vdot, Re, f
    double precision :: dx, dt    
    double precision :: Nu, hc_fb, hc_pf, hc_gw
    double precision :: Data_time_step, Store_dt_min, Store_dt_min_details, Store_dt, Store_dt_details
    double precision :: Total_sim_time,Curr_sim_time, Next_store_time, Next_store_time_details
    double precision :: Ag, Bg, Ap, Bp, Af, Ab, Bb
    double precision :: UsefulHeatTotal, Total_I
    double precision :: Tout_avg_stor_dt, Count_stor_dt
    double precision :: Tg_iter, Tp_iter, Tf_iter, Tb_iter, Col_eff
    double precision :: TOL_delta_Tg, TOL_delta_Tp, TOL_delta_Tf, TOL_delta_Tb
    integer :: Nx, Nt, iterations, iterations_tot, j, output_unit_number, output_DETAILS_unit_number, index, jm1, jp1
    double precision, allocatable :: x_position(:), hr_pg(:), hr_gs(:), hr_pb(:), Ra(:), Nu_nc(:), hnc_pg(:)
    double precision, allocatable :: Tg(:), Tg_new(:), Tp(:), Tp_new(:), Tf(:), Tf_new(:), Tb(:), Tb_new(:)
    double precision, allocatable :: delta_Tg(:), delta_Tp(:), delta_Tf(:), delta_Tb(:) 
    
    character(len=100) :: output_file = 'Output_SAH_Stutt_HOURLYDATA.txt'       ! Name of output file to be used as input by donwstream component
    character(len=100) :: output_file_details = 'Output_SAH_Stutt_HOURLYDATA_details.txt' ! Name of output file with detailed data
    
    output_unit_number = 11
    output_DETAILS_unit_number = 12
    
    pi = 3.14159265359
    sigma = 5.6697e-8           ! Stefan-Bolzmann constant (W/m2-K4)
    relax = 1.0                 ! Underrelaxation factor, can be lower than 1 for stability if required, but more than 0
    Data_time_step = 3600.        ! This is the time resolution of the input data (s)
    
    !Call the subroutine to read data from a start row to an end row
    call read_file("Data_Stuttgart_July-Aug2021_hourly.txt", 2, 120)
    
    ! COLLECTOR PROPERTIES
    L = 6.               ! Collector length (m)
    W = 2.                    ! Collector width (m)
    H = 0.05        ! Gap of air channel (m)
    gap_gp = 0.05        	! Gap between cover and absorber (m) 
    tilt = 25*pi/180     ! Collector tilt (rad) (change only the first number which is the angle in degrees)
    
    thick_g = 0.003; k_g = 1.05; c_g = 660.; rho_g = 2500.0   ! Cover thickness (m), conductivity (W/m-K), specific heat (J/kg-K) and density (kg/m3)
    thick_p = 0.001; k_p = 220.0; c_p = 890.; rho_p = 2700.0   ! Absorber thickness (m), conductivity (W/m-K), specific heat (J/kg-K) and density (kg/m3)
    thick_b = 0.001; k_b = 220.0; c_b = 890.; rho_b = 2700.0   ! Back plate thickness (m), conductivity (W/m-K), specific heat (J/kg-K) and density (kg/m3)
    
    thick_i = 0.1; k_i = 0.04                               ! Back insulation thickness (m) and conductivity (J/kg-K)     
    alpha_g = 0.06; tau_g = 0.84; e_g = 0.90                ! Absorptance, transmittance and emittance of cover (-)
    alpha_p = 0.95; e_p_up = 0.95; e_p_down = 0.95          ! Absorptance and emittance of absorber (-)
    alpha_b = 0.95; e_b = 0.95                              ! Absorptance and emittance of back plate (-)   
    
    ! AIR PROPERTIES, ASSUMED CONSTANT
    rho_f = 1.127                           ! Air density (at about 40  C) (kg/m3)
    c_f = 1007.                             ! Air specific heat (J/kg-K)
    mu_f = 1.907e-5                         ! Air viscosity (at 40  C) (Pa-s)
    k_f = 0.02735                           ! Themal conductivity of air (at 40  C) (W/m-K)
    Pr = 0.706                              ! Prandtl number of air
    
    ! DEFINE MESH SIZE
    dx = 0.1                                    ! Mesh size (tentative) (m)
    dt = 10.                                    ! delta t desired (s)
    
    ! TOLERANCE FOR THE DIFFERENCE BETWEEN CONSECUTIVE ITERATIONS
    TOL_delta_Tg = 0.00001      ! Tolerance of nodal cover temperature ( C)
    TOL_delta_Tp = 0.00001      ! Tolerance of nodal absorber temperature ( C)
    TOL_delta_Tf = 0.00001      ! Tolerance of nodal fluid temperature ( C)
    TOL_delta_Tb = 0.00001      ! Tolerance of nodal bottom plate temperature ( C)
    
    ! CHOOSE HOW OFTEN TO STORE DATA
    Store_dt_min = 1.           ! Data storage interval (min)
    Store_dt_min_details = 5.   ! Data storage interval for detailed collector data (min)
    
    ! START OF CALCULATIONS
    Dh = 4*W*H/(2*(W+H))        ! Hydraulic diameter of flow passage (m)
    Ac = L*W                    ! Collector area (m2)
    Ub = k_i/thick_i            ! Back loss coefficient (W/m2-K)
    th_diff_f = k_f/(c_f*rho_f) ! Air thermal diffusivity (m2/s)
    
    Nx = nint(L/dx)             ! Number of space intervals  
    
    allocate(x_position(0:Nx))
    
    do j=0,Nx
        x_position(j) = j*L/Nx   ! Vector with the position of all nodes, m
    end do
    
    dx = x_position(1)-x_position(0)    ! Actual dx (m)           
    
    ! CALCULATE CONSTANT FACTORS FOR EQUATIONS
    Ag = k_g*thick_g/(dx**2); Bg = rho_g*thick_g*c_g
    Af = 1/(rho_f*c_f*H)
    Ap = k_p*thick_p/(dx**2); Bp = rho_p*thick_p*c_p
    Ab = k_b*thick_b/(dx**2); Bb = rho_b*thick_b*c_b
    
    allocate(hr_pg(0:Nx))           ! Allocating vector of radiation heat transfer coefficients between absorber and cover, which are dependent on Tp and Tg
    allocate(hr_gs(0:Nx))           ! Allocating vector of sky radiation heat transfer coefficients which are dependent on Tg and Ts
    allocate(hr_pb(0:Nx))           ! Allocating vector of radiation heat transfer coefficients betweeen absorber and back plate, which are dependenet on Tp and Tb

    allocate(Ra(0:Nx))              ! Allocating vector of Rayleigh number
    allocate(Nu_nc(0:Nx))           ! Allocating vector of natural convection Nusselt number
    allocate(hnc_pg(0:Nx))          ! Allocating vector of natural convection heat transfer coefficient
    
    Store_dt = Store_dt_min*60  ! Data storage interval, s
    Store_dt_details = Store_dt_min_details*60  ! Data storage interval for detailed collector data, s
    
    Total_sim_time = size(time_row)*Data_time_step      ! The total simulation time, s
      
    ! Initialize main temperature vectors
    
    allocate(Tg(0:Nx)); Tg = Tamb(0)            ! Defines the size of the vector Tg and initializes to Tamb(0)
    allocate(Tg_new(0:Nx)); Tg_new = Tamb(0)    ! Defines the size of the vector Tg_new and initializes to Tamb(0)
 
    allocate(Tp(0:Nx)); Tp = Tamb(0) + 0.1               ! Defines the size of the vector Tp and initializes to Tamb(0)+0.1
    allocate(Tp_new(0:Nx)); Tp_new = Tamb(0) + 0.1       ! Defines the size of the vector Tp_new and initializes to Tamb(0)+0.1
    
    allocate(Tf(0:Nx)); Tf = Tamb(0)                 ! Defines the size of the vector Tf and initializes to Tamb(0)
    allocate(Tf_new(0:Nx)); Tf_new = Tamb(0)       ! Defines the size of the vector Tf_new and initializes to Tamb(0)
 
    allocate(Tb(0:Nx)); Tb = Tamb(0)               ! Defines the size of the vector Tb and initializes to Tamb(0)
    allocate(Tb_new(0:Nx)); Tb_new = Tamb(0)       ! Defines the size of the vector Tb_new and initializes to Tamb(0)
  
    UsefulHeatTotal = 0.            ! variable to keep track of useful heat over ambient T
    Total_I = 0.                    ! variable to keep track of total solar radiation (J)

    Curr_sim_time = 0.0         ! The current simulation time during the calculations, s
    Next_store_time = 0.                ! Initializing the variable keeping track of the next simulation time to store values
    Next_store_time_details = 0.        ! Initializing the variable keeping track of the next simulation time to store detail values
    
    ! HERE WE WRITE TO THE FILE THAT IS GOING TO BE USED AS INPUT BY THE DOWNSTREAM COMPONENT
    open(unit=output_unit_number, file=trim(output_file), status='replace', action='write')
    
    write(output_unit_number, '(11(A,F8.3))') "L:", L, " W:", W, " H:", H, " gap_gp:", gap_gp, " tilt:", tilt, &
        " thick_g:", thick_g, " thick_p:", thick_p, " thick_b:", thick_b, " k_g:", k_g, " k_p:", k_p, " k_b:", k_b
    write(output_unit_number, '(6(A,F8.1),2(A,F8.3))') "c_g:", c_g, " c_p:", c_p, " c_b:", c_b, " rho_g:", rho_g, " rho_p:", rho_p, &
        " rho_b:", rho_b, " thick_i:", thick_i, " k_i:", k_i
    write(output_unit_number, '(8(A,F8.3))') "alpha_g:", alpha_g, " alpha_p:", alpha_p, " alpha_b:", alpha_b, " tau_g:", tau_g, " e_g:", e_g, &
        " e_p_up:", e_p_up, " e_p_down:", e_p_down, " e_b:", e_b
    write(output_unit_number, '(A,F8.2,A,F8.1)') "dx:", dx, " dt:", dt
    write(output_unit_number, '(A,A12,F9.3,F8.1,F8.1,F8.1,F8.1,F8.0,F8.4,F13.3,F11.8)') Date(0), Hour_min(0), Hour_dec(0), Curr_sim_time/60., &
        Ins(0), Tamb(0), wind(0), Patm(0), Xamb(0), Tf(Nx), mdot_a(0)
    write(output_unit_number, '(A,A12,A9,A8,A8,A8,A8,A8,A8,A13,A11)') "Date", "Hour_min", "Hour_dec", "t_min", "Ins", "Tamb", "Wind", "Patm", "X", &
        "Tpcm_out_avg", "mdot_a"
    
    ! HERE WE WRITE TO THE FILE THAT IS GOING TO CONTAIN DETAILED TEMPERATURE IN SAH
    open(unit=output_DETAILS_unit_number, file=trim(output_file_details), status='replace', action='write')
    
    write(output_DETAILS_unit_number, '(11(A,F8.3))') "L:", L, " W:", W, " H:", H, " gap_gp:", gap_gp, " tilt:", tilt, &
        " thick_g:", thick_g, " thick_p:", thick_p, " thick_b:", thick_b, " k_g:", k_g, " k_p:", k_p, " k_b:", k_b
    write(output_DETAILS_unit_number, '(6(A,F8.1),2(A,F8.3))') "c_g:", c_g, " c_p:", c_p, " c_b:", c_b, " rho_g:", rho_g, " rho_p:", rho_p, &
        " rho_b:", rho_b, " thick_i:", thick_i, " k_i:", k_i
    write(output_DETAILS_unit_number, '(8(A,F8.3))') "alpha_g:", alpha_g, " alpha_p:", alpha_p, " alpha_b:", alpha_b, " tau_g:", tau_g, " e_g:", e_g, &
        " e_p_up:", e_p_up, " e_p_down:", e_p_down, " e_b:", e_b
    write(output_DETAILS_unit_number, '(A,F8.2,A,F8.1)') "dx:", dx, " dt:", dt
    write(output_DETAILS_unit_number, '(A,A12,A9,A8,A8,A8,A8)', advance = 'no') "Date", "Hour_min", "Hour_dec", "t_min", "Ins", "Tamb", "Wind"
    
    do j = 0,size(x_position)-1
        write(output_DETAILS_unit_number, '(F8.2)', advance = 'no') x_position(j)
    end do
    do j = 0,size(x_position)-1
        write(output_DETAILS_unit_number, '(F8.2)', advance = 'no') x_position(j)
    end do
    do j = 0,size(x_position)-1
        write(output_DETAILS_unit_number, '(F8.2)', advance = 'no') x_position(j)
    end do
    do j = 0,size(x_position)-1
        write(output_DETAILS_unit_number, '(F8.2)', advance = 'no') x_position(j)
    end do
    
    write(output_DETAILS_unit_number, '(A)') " Eff"
    
    write(output_DETAILS_unit_number, '(A,A12,F9.3,F8.1,F8.1,F8.1,F8.1)', advance = 'no') Date(0), Hour_min(0), Hour_dec(0), Curr_sim_time/60., &
        Ins(0), Tamb(0), wind(0)
    
    do j = 0,size(x_position)-1
        write(output_DETAILS_unit_number, '(F8.2)', advance = 'no') Tg(j)
    end do
    do j = 0,size(x_position)-1
        write(output_DETAILS_unit_number, '(F8.2)', advance = 'no') Tp(j)
    end do
    do j = 0,size(x_position)-1
        write(output_DETAILS_unit_number, '(F8.2)', advance = 'no') Tf(j)
    end do   
    do j = 0,size(x_position)-1
        write(output_DETAILS_unit_number, '(F8.2)', advance = 'no') Tb(j)
    end do
    write(output_DETAILS_unit_number,*)
    
    Next_store_time = Next_store_time + Store_dt            ! set next storage time
    Next_store_time_details = Next_store_time_details + Store_dt_details  ! set next storage time in detailed output file

    Tout_avg_stor_dt = 0.0
    Count_stor_dt = 1              ! This counts how many timesteps have passed before storing results in the first output file, so that 
                                   ! the avg of outlet Ta during the store_dt can be calculated. Starts at 1 because the least is to write data every time step.

    allocate(delta_Tg(0:Nx)); delta_Tg = 1.       ! Matrix of iteration changes of Tg in each node, inititalized to a value above the tolerance desired
    allocate(delta_Tp(0:Nx)); delta_Tp = 1.       ! Matrix of iteration changes of Tp in each node, initialized to a value above the tolerance desired
    allocate(delta_Tf(0:Nx)); delta_Tf = 1.       ! Matrix of iteration changes of Tf in each node, inititalized to a value above the tolerance desired
    allocate(delta_Tb(0:Nx)); delta_Tb = 1.       ! Matrix of iteration changes of Tb in each node, initialized to a value above the tolerance desired

    iterations_tot = 0    
    
    ! TIME LOOP STARTS HERE    
    do while(Curr_sim_time < Total_sim_time)
        iterations = 0
        ! Retrieve needed input data from weather file                        
        h_dec = Hour_dec(int(Curr_sim_time/Data_time_step))
        I = Ins(int(Curr_sim_time/Data_time_step))
        Ta = Tamb(int(Curr_sim_time/Data_time_step))
        v_wind = wind(int(Curr_sim_time/Data_time_step))
        p_atm = Patm(int(Curr_sim_time/Data_time_step))
        Xa = Xamb(int(Curr_sim_time/Data_time_step))
        Tdp = Tdewp(int(Curr_sim_time/Data_time_step))
        mdot_f = mdot_a(int(Curr_sim_time/Data_time_step))
            
        Vdot = mdot_f/rho_f                 ! Airflow rate (m3/s)
        v = Vdot/(W*H)                      ! Average air speed in air channel (m/s)
        spec_Vdot = Vdot/Ac                 ! Airflow rate per unit collector area (m3/s-m2)
        Re = v*Dh*rho_f/mu_f                ! Reynolds number for flow in collector
        
    
        if (Re <= 2550.0) then              ! Calculates the Nusselt number based on flow conditions
            Nu = 5.385 + 0.148*Re*H/L
        else if (Re <= 10000.0) then
            Nu = 4.4e-4*Re**1.2 + 9.37*Re**0.471*H/L
        else
            Nu = 0.03*Re**0.74 + 0.788*Re**0.74*H/L
        end if
    
        hc_fb = Nu*k_f/Dh               ! Conv. heat transf. coeff. fluid-back plate (W/m2-K)
        hc_pf = Nu*k_f/Dh               ! Conv. heat transf. coeff. absorber-fluid (W/m2-K)
        
        Ts = ((Ta + 273.15)*(0.711 + 0.0056*Tdp + 0.000073*Tdp**2 + 0.013*cos(15*h_dec*pi/180))**0.25) - 273.15 ! Sky temperature ( C) (Duffie and Beckmann)
                                                                                    ! This equation makes use of the hour of the day (h_dec in the equation)
                                                                                    ! Here the h_dec is the hour of day in decimal, for example 6:30 am would be 6.5
                                                                                    ! and 2:15 pm would be 14.25
        
        hc_gw = 5.7 + 3.8*v_wind        ! Wind conv. heat trans. coeff. (W/m2-K)
        
        hr_pg = sigma*((Tp+273)**2 + (Tg+273)**2)*((Tp+273) + (Tg+273))/(1/e_p_up + 1/e_g - 1)      ! Radiation heat transfer coefficients (W/m2-K)
        hr_pb = sigma*((Tp+273)**2 + (Tb+273)**2)*((Tp+273) + (Tb+273))/(1/e_p_down + 1/e_b - 1)
        hr_gs = sigma*e_g*((Tg+273)**2 + (Ts+273)**2)*((Tg+273) + (Ts+273))
        
        Ra = 9.81*(1/(Tf + 273))*abs(Tp - Tg)*(gap_gp**3)*rho_f/(mu_f*th_diff_f)         ! Rayleigh number
        Nu_nc = 1 + 1.44*(1 - 1708*sin(1.8*tilt)**1.6/(Ra*cos(tilt)))*max(1 - 1708/(Ra*cos(tilt)), 0.0) + max((Ra*cos(tilt)/5830)**(1./3.) - 1, 0.0) ! Nussel number for natural convection (Duffie and Beckman)
        
        hnc_pg = Nu_nc*k_f/gap_gp
        
        do while (maxval(delta_Tg) > TOL_delta_Tg .or. maxval(delta_Tp) > TOL_delta_Tp .or. maxval(delta_Tf) > TOL_delta_Tf .or. maxval(delta_Tb) > TOL_delta_Tb) !# This is the iteration control, the tolerance value can be varied
            
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
                
                Tg_iter = (Bg*Tg(j) + dt*(Ag*(Tg_new(jm1) + Tg_new(jp1)) + (hnc_pg(j) + hr_pg(j))*Tp_new(j) + alpha_g*I + hc_gw*Ta + hr_gs(j)*Ts))/(Bg + dt*(2*Ag + hnc_pg(j) + hc_gw + hr_pg(j) + hr_gs(j)))
            
                Tg_iter = relax*Tg_iter + (1-relax)*Tg_new(j)             ! Applying underrelaxation 
            
                delta_Tg(j) = abs(Tg_iter - Tg_new(j))          ! Recording the change in temperature from past iteration
                        
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
            
                Tp_iter = (Bp*Tp(j) + dt*((hnc_pg(j)+ hr_pg(j))*Tg_new(j) + Ap*(Tp_new(jm1) + Tp_new(jp1)) + hc_pf*Tf_new(j) + hr_pb(j)*Tb_new(j) + alpha_p*tau_g*I))/(Bp + dt*(2*Ap + hnc_pg(j) + hc_pf + hr_pg(j) + hr_pb(j)))
    
                Tp_iter = relax*Tp_iter + (1-relax)*Tp_new(j)             ! Applying underrelaxation 
            
                delta_Tp(j) = abs(Tp_iter - Tp_new(j))          ! Recording the change in temperature from past iteration
            
                Tp_new(j) = Tp_iter
            end do
            
            ! LOOP OVER AIR NODES FROM i=0 TO Nx
            Tf_iter = Ta
            Tf_new(0) = Tf_iter
            delta_Tf(0) = Tf_iter - Tf_new(0)
            do j = 1,Nx
                Tf_iter = (Tf(j) + dt*(Af*hc_pf*Tp_new(j) + (v/dx)*Tf_new(j-1) + Af*hc_fb*Tb_new(j)))/(1 + dt*(v/dx + Af*(hc_pf + hc_fb)))
            
                Tf_iter = relax*Tf_iter + (1-relax)*Tf_new(j)             ! Applying underrelaxation 
            
                delta_Tf(j) = abs(Tf_iter - Tf_new(j))          ! Recording the change in temperature from past iteration
            
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
            
                delta_Tb(j) = abs(Tb_iter - Tb_new(j))          ! Recording the change in temperature from past iteration
            
                Tb_new(j) = Tb_iter
            end do
        
            iterations = iterations + 1
            iterations_tot = iterations_tot + 1
           
        end do            
            
        print *, "max_delta_Tg", maxval(delta_Tg)
        print *, "max_delta_Tp", maxval(delta_Tp)
        print *, "max_delta_Tf", maxval(delta_Tf)
        print *, "max_delta_Tb", maxval(delta_Tb)
        
        print *, "t", Curr_sim_time
        print *, "iterations", iterations
        
        delta_Tg = 1.0      ! Reset the values of the deltas in each node to a value above the tolerance desired
        delta_Tp = 1.0
        delta_Tf = 1.0
        delta_Tb = 1.0
        
        Tg = Tg_new         ! Replacing the new calculated temperatures into the vectors of current temperature to move into next timestep
        Tp = Tp_new
        Tf = Tf_new
        Tb = Tb_new
        
        if (Tf(Nx) > Ta) then
            UsefulHeatTotal = UsefulHeatTotal + (Tf(Nx)-Ta)*mdot_f*c_f*dt
        end if
    
        Total_I = Total_I + I*Ac*dt
        
        Tout_avg_stor_dt =  (Tout_avg_stor_dt*(Count_stor_dt - 1) + Tf(Nx))/Count_stor_dt    ! This is used to calculate the average outlet T during the store_dt.
                                                                                         ! It calculates the mean of all Tout values between store_dts to give a value
                                                                                         ! of mean Tout for the whole storage period. For example if i store data every minute
                                                                                         ! and my ts is 10s, the code calculates a Tout every 10s, so here we calculate the mean
                                                                                         ! of 6 Tout values which represent the mean Tout of that minute.                      
    
        Curr_sim_time = Curr_sim_time + dt      ! Advance time step
        
        
        ! STORING DATA IN OUTPUT FILES
        index = int(Curr_sim_time/Data_time_step)

        if (mod(Curr_sim_time, Data_time_step) == 0) then
            index = index - 1
        end if
    
        if (Next_store_time - Curr_sim_time < dt) then
            write(output_unit_number, '(A,A12,F9.3,F8.1,F8.1,F8.1,F8.1,F8.0,F8.4,F13.3,F11.8)') Date(index), Hour_min(index), Hour_dec(index), &
                Curr_sim_time/60., Ins(index), Tamb(index), wind(index), Patm(index), Xamb(index), Tout_avg_stor_dt, mdot_a(index)
            print *, "saving..."
            
            Tout_avg_stor_dt = 0
            Count_stor_dt = 0
            Next_store_time = Next_store_time + Store_dt
        end if
    
        if (Next_store_time_details - Curr_sim_time < dt) then
            write(output_DETAILS_unit_number, '(A,A12,F9.3,F8.1,F8.1,F8.1,F8.1)', advance = 'no') Date(index), Hour_min(index), Hour_dec(index), Curr_sim_time/60., &
            Ins(index), Tamb(index), wind(index)
    
            do j = 0,size(x_position)-1
            write(output_DETAILS_unit_number, '(F8.2)', advance = 'no') Tg(j)
            end do
            do j = 0,size(x_position)-1
                write(output_DETAILS_unit_number, '(F8.2)', advance = 'no') Tp(j)
            end do
            do j = 0,size(x_position)-1
                write(output_DETAILS_unit_number, '(F8.2)', advance = 'no') Tf(j)
            end do   
            do j = 0,size(x_position)-1
                write(output_DETAILS_unit_number, '(F8.2)', advance = 'no') Tb(j)
            end do
            if (I > 0) then
                write(output_DETAILS_unit_number, '(F8.3)') ((Tf(Nx)-Ta)*(mdot_f)*c_f/(I*Ac))
            else
                write(output_DETAILS_unit_number,*)
            end if
            Next_store_time_details = Next_store_time_details + Store_dt_details
        end if
    
        Count_stor_dt = Count_stor_dt + 1
    end do
    
    Col_eff = UsefulHeatTotal/Total_I
    write(output_DETAILS_unit_number, '(A,F15.2)') "Total I (J):", Total_I
    write(output_DETAILS_unit_number, '(A,F15.2)') "Total Useful heat (J):", UsefulHeatTotal 
    
    write(output_DETAILS_unit_number, '(A,F15.3)') "Collector efficiency:", Col_eff
    
    deallocate(Date,Hour_min,Hour_dec,time_row,Ins,Tamb,rh_amb,wind,Patm,Xamb,Tdewp,mdot_a,x_position)
    deallocate(Tg,Tg_new,Tp,Tp_new,Tf,Tf_new,Tb,Tb_new,hr_pg,hr_gs,hr_pb,Ra,Nu_nc,hnc_pg)
    deallocate(delta_Tg,delta_Tp,delta_Tf,delta_Tb)

    end program SAH