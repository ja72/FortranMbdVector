!****************************************************************************
!
!  PROGRAM: FortranMbdVector
!
!  PURPOSE: Rigid Body Simulation using custom vector types
!
!****************************************************************************

    program FortranMbdVector
    use mod_vectors_simulation
    implicit none
    
    integer, parameter :: report_count = 36
    
    ! Real parameters
    type(vector3), parameter :: &
        gravity = vector3(0.0_wp, 0.0_wp, -10.0_wp)      ! gravity [m/s^2]

    ! Variables
    integer :: lu, du, ios
    
    !call mbd_solver(9*36)
    call mbd_solver(360*2250, 2.0_wp)
    
    contains
    
    subroutine mbd_solver(n_steps, end_time, angularResolutionDegrees)
    integer, intent(in) :: n_steps
    real(wp), intent(in), optional :: end_time, angularResolutionDegrees
    
    ! Cylinder: Length=90 [mm], Radius=30 [mm]
    real(wp), parameter :: L = 90*mm
    real(wp), parameter :: D = 60*mm
    real(wp), parameter :: &        
        eps = 0.87_wp, mu=0.15_wp                       ! COR and COF           
    
    type(vector3) :: pos, vee, omg
    type(quaternion) :: ori
    real(wp), allocatable :: Y(:,:)
    real(wp) :: t_end, h, theta
    type(contact_plane) :: floor
    type(simulation) :: mbd
    logical :: itsopen
    
    integer(8) :: tic, toc, cpurate
    real(real32) :: cpu, ops
    character(1) :: key
    integer, parameter :: ESC = 27

    class(rigid_body), allocatable :: body
        
    print *, "Fortran", n_steps, "steps."
    print *, ""
        
    body = rigid_body( shapes(cylinder, [L, D/2]), mass=2.0_wp)
    
    pos = (L/2)*k_                              ! initial position [m]
    ori = q_eye                                 ! initial orientation (quaternion)
    omg = 5.0_wp*j_                             ! initial rotational velocity [rad/s]
    vee = (-1.0_wp)*k_                          ! initial velocity [m/s]
    
    call body%set_pose(pos, ori)
    call body%set_motion(vee, omg)
    
    floor = contact_plane(o_, k_, eps, mu)
    
    mbd = simulation(body,floor,gravity,sim_step_taken)
    
    if(present(end_time)) then
        t_end = end_time
    else
        t_end = 2.0_wp
    end if
    if(present(angularResolutionDegrees)) then
        theta = angularResolutionDegrees
    else
        theta = 0.5_wp
    end if
    
    h = body%est_max_time_step(n_steps, t_end, theta)
        
    t_end = n_steps * h    
    
    if(.false.) then
        open(newunit=lu, file='results-f.csv', status='replace', form='formatted', iostat=ios)
    else
        lu = 9999
    end if
    
    inquire(unit=lu, opened=itsopen)     
    
    if(itsopen .and. ios/=0) then
        error stop 'File Creation Error'
    end if
            
    call write_table_header()
    if(itsopen) then
        call write_csv_header(lu)
    end if
        
    call SYSTEM_CLOCK(tic,cpurate)    
    call mbd%run(t_end, n_steps)    
    call SYSTEM_CLOCK(toc,cpurate)
    
    if(itsopen) then
        close(lu)
    end if
            
    cpu = real(toc-tic)/cpurate
    ops = n_steps / cpu
    
    print *, ""
    print '(1x,"steps=",g0,1x,"time=",g0,1x,"kops=",g0)', n_steps, cpu, ops/1000

    end subroutine
    
    subroutine sim_step_taken(sim,n)
    class(simulation), intent(in) :: sim
    integer, intent(in) :: n
    integer :: i
    real(wp) :: t, Y(state_size)
    logical itsopen 
    type(contact_plane) :: ctx
    
        i =sim%step
        if( i==n .or. mod(i, n/report_count) ==0 ) then
            call write_table_row(sim)
        end if    
        
        if(itsopen) then
            call write_csv_row(sim,lu)
        end if
    end subroutine
    
    subroutine write_table_header(io)
    integer, intent(in), optional :: io
    integer :: unit
    character(len=*), parameter :: fmt = '(a9,1x,a9,"|",3a9,"|",4a9,"|",3a9,"|",3a10,"|",a8)'
    
    if( present(io) ) then
        unit = io
    else
        unit = 6
    end if
    
    write (unit, fmt) &
        "#", "time", &
        "pos-x", "pos-y","pos-z", &
        "ori-x", "ori-y", "ori-z", "ori-w", &
        "mom-x", "mom-y", "mom-z", &
        "ang-x", "ang-y", "ang-z", "J"

    write (unit, fmt) &
        "", "[s]", &
        "[m]", "[m]","[m]", &
        "[]", "[]", "[]", "[]", &
        "[kg m/s]", "[kg m/s]", "[kg m/s]", &
        "[g m^2/s]", "[g m^2/s]", "[g m^2/s]", &
        "[N s]"
    
    write (unit, fmt) &
        "---", "---", &
        "---", "---","---", &
        "---", "---", "---", "---", &
        "---", "---", "---", &
        "---", "---", "---", "---"
    end subroutine
    
    subroutine write_table_row(sim, io)
    type(simulation), intent(in) :: sim
    integer, intent(in), optional :: io
    real(wp) :: t, Y(state_size), Jn
    integer :: unit, i
    character(len=*), parameter :: fmt = '(i9,1x,f9.3,"|", 3f9.5,"|", 4f9.4,"|",3f9.5,"|", 3f10.4,"|", f8.4)'        
    
    if( present(io) ) then
        unit = io
    else
        unit = 6
    end if
    
    i = sim%step
    t = sim%time
    Y = sim%current
    Jn = sim%contact%Jn
    
    write (unit, fmt) i, t, Y(1:3), Y(4:7), Y(8:10), 1000*Y(11:13), Jn
    end subroutine

    subroutine write_csv_header(io)
    integer, intent(in) :: io
    integer :: unit
    character(len=*), parameter :: fmt = '(g0,*(",",g0))'
    
    unit = io
    
    write (unit, fmt) &
        "#", "time [s]", &
        "pos-x [m]", "pos-y [m]","pos-z [m]", &
        "ori-x []", "ori-y []", "ori-z []", "ori-w []", &
        "mom-x [kg m/s]", "mom-y [kg m/s]", "mom-z [kg m/s]", &
        "ang-x [g m^2/s]", "ang-y [g m^2/s]", "ang-z [g m^2/s]", &
        "J [N s]"
    
    end subroutine
    
    subroutine write_csv_row(sim, io)
    type(simulation), intent(in) :: sim
    integer, intent(in) :: io
    real(wp) :: t, Y(state_size), Jn
    integer :: unit, i
    character(len=*), parameter :: fmt = '(g0,*(",",g0))'
    
    unit = io
    
    i = sim%step
    t = sim%time
    Y = sim%current
    Jn = sim%contact%Jn
    
    write (unit, fmt) i, t, Y(1:3), Y(4:7), Y(8:10), 1000*Y(11:13), Jn
    end subroutine

    end program FortranMbdVector

