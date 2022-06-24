    module mod_vectors_dynamics
    use mod_vectors_quaternion
    use mod_vectors_geometry
    implicit none
    
    ! Real parameters
    real(wp), parameter :: &        
        steel_density = 7860.0_wp, &                    ! shape density [kg/m^3]
        eps = 0.87_wp, mu=0.15_wp                       ! COR and COF           
    
    
    type(vector3), parameter :: &
        gravity = vector3(0.0_wp, 0.0_wp, -10.0_wp)      ! gravity [m/s�]
    
    ! Simulation parameters
    integer, parameter :: state_size = 13, report_count = 36
        
    ! Rigid Body Definition
    type :: rigid_body
        class(shapes), allocatable :: shape
        real(wp) :: mass, mmoi(3)
        type(vector3) :: cg                         ! cg from ref. point
        type(vector3) :: initial_position           ! at ref. point
        type(quaternion) :: initial_orientation             
        type(vector3) :: initial_velocity           ! at ref. point
        type(vector3) :: initial_omega    
    contains
        procedure :: initial_state => rb_get_initial_state
        procedure :: set_pose => rb_set_cg_pose
        procedure :: set_motion => rb_set_cg_motion
        procedure :: get_motion => rb_get_motion
        procedure :: integrate => rb_integrate
        procedure :: simulate => rb_simulate
    end type
    
    interface rigid_body
        module procedure :: rb_from_shape, rb_new
    end interface
    
    type :: contact_plane
        real(wp) :: epsilon, mu
        type(vector3) :: origin
        type(vector3) :: normal
        type(vector3) :: slip
        real(wp) :: delta, v_imp, Jn
        real(wp) :: v_slip, Js
        logical :: active
    contains
        procedure :: contact => ctx_contact_calc    
        procedure :: reset => ctx_reset
    end type
        
    interface contact_plane
        module procedure :: ctx_new_plane
    end interface
    
    interface 
        subroutine step_handler(i,n,t,Y,ctx)
        import
        integer, intent(in) :: i,n
        real(wp), intent(in) :: t, Y(:)
        type(contact_plane), intent(in) :: ctx
        end subroutine
    end interface
    
    character(len=*), parameter :: fmth = '(a9  ,"|", 3a9  ,"|", 4a9  ,"|",3a9  ,"|", 3a10  ,"|", a8)'
    character(len=*), parameter :: fmtv = '(f9.3,"|", 3f9.5,"|", 4f9.4,"|",3f9.5,"|", 3f10.6,"|", f8.4)'    
        
    contains
    
    function rb_new(s,mass,cg,mmoi) result(rb)
    class(shapes), intent(in) :: s
    type(rigid_body) :: rb
    real(wp), intent(in) :: mass, mmoi(3)
    type(vector3), intent(in) :: cg
        rb%shape = s
        rb%mass = mass
        rb%cg = cg
        rb%mmoi = mmoi
        rb%initial_position = o_
        rb%initial_orientation = q_eye
        rb%initial_velocity = o_
        rb%initial_omega = o_
    end function
    
    function rb_from_shape(s, density) result(rb)
    class(shapes), intent(in) :: s
    real(wp), intent(in), optional :: density
    type(rigid_body) :: rb

    if( present(density) ) then
        rb = rigid_body(s, density * s%volume, s%center, density * s%vmmoi)
    else
        rb = rigid_body(s, 1000.0_wp * s%volume, s%center, 1000.0_wp * s%vmmoi)
    end if

    end function
    
    
    pure function rb_get_initial_state(rb) result(Y)
    class(rigid_body), intent(in) :: rb
    real(wp):: Y(state_size)
        Y = rb_get_state(rb, &
            rb%initial_position, &
            rb%initial_orientation, &
            rb%initial_velocity, &
            rb%initial_omega)
    end function
    
    pure subroutine rb_set_cg_pose(rb, pos, ori)
    class(rigid_body), intent(inout) :: rb
    type(vector3), intent(in) :: pos
    type(quaternion), intent(in) :: ori
    type(vector3) :: c
    type(matrix3) :: R
        R = rot(ori)
        c = R*rb%cg
        rb%initial_position = pos - c
        rb%initial_orientation = ori
    end subroutine
    
    pure subroutine rb_set_cg_motion(rb, vee, omg)
    class(rigid_body), intent(inout) :: rb
    type(vector3), intent(in) :: vee, omg
    type(vector3) :: c
    type(matrix3) :: R
        R = rot(rb%initial_orientation)
        c = R*rb%cg
        rb%initial_velocity = vee + cross(c,omg)
        rb%initial_omega = omg
    end subroutine
    
    pure function rb_get_state(rb, pos_b, ori, vee_b, omg) result(Y)
    class(rigid_body), intent(in) :: rb
    type(vector3), intent(in) :: pos_b, vee_b, omg
    type(quaternion), intent(in) :: ori
    real(wp):: Y(state_size)
    type(vector3) :: p, L_b, c
    type(matrix3) :: R, I
        R = rot(ori)
        c = R * rb%cg
        I = R*diag(rb%mmoi)*R%transpose()
        p = rb%mass*(vee_b - cross(c, omg))
        L_b = I*omg + cross(c, p)        
        Y = [ &
            pos_b%x, pos_b%y, pos_b%z, &
            ori%vector%x, ori%vector%y, ori%vector%z, ori%scalar, &
            p%x, p%y, p%z, &
            L_b%x, L_b%y, L_b%z]
    end function
    
    pure subroutine rb_get_motion(rb, Y, vee_b, omg)
    class(rigid_body), intent(in) :: rb
    real(wp), intent(in) :: Y(state_size)
    type(vector3), intent(out) :: vee_b
    type(vector3), intent(out) :: omg
    type(vector3) :: c, p, L_b
    type(quaternion) :: ori
    type(matrix3) :: R, I_inv
    
        ori = Y(4:7)
        p = Y(8:10)
        L_b = Y(11:13)
        R = rot(ori)
        c = R * rb%cg
        I_inv = R*diag(1/rb%mmoi)*R%transpose()
        omg = I_inv*(L_b - cross(c, p))
        vee_b = p/rb%mass - cross(omg, c)
    end subroutine
    
    pure subroutine y_normalize(Y)
    real(wp), intent(inout) :: Y(state_size)
        ! Normalize the quaternion in the state vector
        Y(4:7) = Y(4:7)/norm2(Y(4:7))
    end subroutine
    
    pure function rb_get_rate(rb, t, Y) result(Yp)
    class(rigid_body), intent(in) :: rb
    real(wp), intent(in) :: t, Y(state_size)
    real(wp) :: Yp(state_size)
    type(vector3) :: F, tau_b, omg, vee_b, pos_b
    type(quaternion) :: ori, qp
    type(vector3) :: c, p, L_b, tau
    type(matrix3) :: R, I_inv
    
        call y_normalize(Y)
    
        pos_b = Y(1:3)
        ori = Y(4:7) 
        p = Y(8:10)
        L_b = Y(11:13)
        R = rot(ori)
        c = R * rb%cg
        I_inv = R*diag(1/rb%mmoi)*R%transpose()
        omg = I_inv*(L_b - cross(c, p))
        vee_b = p/rb%mass - cross(omg, c)
        qp = 0.5_wp*(quaternion(omg) .o. ori)
        F = rb%mass * gravity
        tau_b = cross(c, F)
        tau = tau_b - cross(vee_b, p)
        Yp = [ &
            vee_b%x, vee_b%y, vee_b%z, &
            qp%vector%x, qp%vector%y, qp%vector%z, qp%scalar, &
            F%x, F%y, F%z, &
            tau%x, tau%y, tau%z]
    
    end function
    
    pure function rb_integrate(rb, t, Y, h) result(Y_next)
    class(rigid_body), intent(in) :: rb
    real(wp), intent(in) :: t, Y(state_size), h
    real(wp) :: Y_next(state_size)
    real(wp) :: k0(state_size),k1(state_size),k2(state_size),k3(state_size)
    
        k0 = rb_get_rate(rb, t, Y)
        k1 = rb_get_rate(rb, t+h/2, Y + (h/2)*k0)
        k2 = rb_get_rate(rb, t+h/2, Y + (h/2)*k1)
        k3 = rb_get_rate(rb, t+h, Y + (h)*k2)
        
        Y_next = Y + (h/6)*(k0+2*k1+2*k2+k3)
        call y_normalize(Y_next)
    end function
    
    pure function ctx_new_plane(origin,normal,epsilon,mu) result(ctx)
    type(contact_plane) :: ctx
    real(wp),intent(in) :: epsilon, mu
    type(vector3),intent(in) :: origin
    type(vector3),intent(in) :: normal    
        ctx%epsilon = epsilon
        ctx%mu = mu
        ctx%origin = origin
        ctx%normal = normal
        ctx%slip   = o_
        ctx%delta  = 0.0_wp
        ctx%v_imp  = 0.0_wp
        ctx%Jn     = 0.0_wp
        ctx%v_slip = 0.0_wp
        ctx%Js     = 0.0_wp
        ctx%active = .false.
    end function
    
    pure subroutine ctx_reset(ctx)
    class(contact_plane), intent(inout) :: ctx
        ctx%slip   = o_
        ctx%delta  = 0.0_wp
        ctx%v_imp  = 0.0_wp
        ctx%Jn     = 0.0_wp
        ctx%v_slip = 0.0_wp
        ctx%Js     = 0.0_wp
        ctx%active = .false.
    end subroutine
        
    function ctx_contact_calc(ctx, rb, Y) result(Y_next)
    class(contact_plane), intent(inout) :: ctx
    type(rigid_body), intent(in) :: rb
    real(wp), intent(in) :: Y(state_size)
    real(wp) :: Y_next(state_size)
    type(vector3) :: p, L_b, omg, vee_b, pos_b
    type(quaternion) :: ori, qp
    type(vector3) :: c, r_A, c_A, v_A, J_imp
    type(matrix3) :: R, I_inv, cx, M_inv
    real(wp) :: s, e, m_eff, m_slip
    logical :: active
            
        ! Body motions
        pos_b = Y(1:3)
        ori = Y(4:7) 
        p = Y(8:10)
        L_b = Y(11:13)
        R = rot(ori)
        c = R * rb%cg
        I_inv = R*diag(1/rb%mmoi)*R%transpose()
        omg = I_inv*(L_b - cross(c, p))
        vee_b = p/rb%mass - cross(omg, c)
        
        ! contact point on surface of shape
        r_A = pos_b + c + rb%shape%nearest_point(ctx%normal, R)
        ctx%delta = dot(ctx%normal, r_A - ctx%origin)
        ! clip contact point on plane
        r_A = r_A - min(ctx%delta,0.0_wp)*ctx%normal
                
        ! get contact distance and speed
        v_A = vee_b + cross(omg, r_A - pos_b)
        ctx%v_imp = dot(ctx%normal, v_A)
        
        ctx%active = ctx%delta<=0 .and. ctx%v_imp < 0
        
        if(ctx%active) then
            ! location of CG relative to contact point
            c_A = pos_b + c - r_A
            cx = cross(c_A)
            ! move inv inertia to contact point
            M_inv = (1/rb%mass) * eye_ - cx*I_inv*cx
            m_eff = 1/dot(ctx%normal, M_inv*ctx%normal)
            ctx%Jn = -(1+ctx%epsilon)*m_eff*ctx%v_imp
            
            ctx%slip = v_A - ctx%v_imp * ctx%normal
            ctx%v_slip = norm2(ctx%slip)
            if( abs(ctx%v_slip) < tiny ) then
                ctx%slip = o_
                ctx%v_slip = 0.0_wp
                ctx%Js = 0.0_wp
            else
                ctx%slip = -ctx%slip/ctx%v_slip
                m_slip = 1/dot(ctx%slip, M_inv*ctx%slip)
                ctx%Js = m_slip*ctx%v_slip
                ctx%Js = min(ctx%mu*ctx%Jn, ctx%Js)
            end if
                        
            J_imp = ctx%Jn*ctx%normal + ctx%Js*ctx%slip
            p = p + J_imp
            L_b = L_b + cross(r_A-pos_b, J_imp)
            
            Y_next = [ &
                pos_b%x, pos_b%y, pos_b%z, &
                ori%vector%x, ori%vector%y, ori%vector%z, ori%scalar, &
                p%x, p%y, p%z, &
                L_b%x, L_b%y, L_b%z]
        else
            ctx%slip   = o_
            ctx%Jn     = 0.0_wp
            ctx%v_slip = 0.0_wp
            ctx%Js     = 0.0_wp
            Y_next = Y
        end if
    end function        
    
    function rb_simulate(rb, t_end, n_steps, update) result(results)
    class(rigid_body), intent(in) :: rb
    integer, intent(in) :: n_steps
    real(wp), intent(in) :: t_end
    real(wp), allocatable :: results(:, :)
    procedure(step_handler), intent(in), pointer, optional :: update
    real(wp) :: h, t, h_next
    real(wp) :: Y(state_size), Y_next(state_size)    
    type(contact_plane) :: floor
    integer :: step
    logical :: report
        
        report = present(update)
        
        allocate(results(0:n_steps, state_size+2))
        
        floor = contact_plane(o_, k_, eps, mu)
        h = t_end/n_steps
        t = 0.0_wp
        Y_next = rb_get_initial_state(rb)
        Y = floor%contact(rb, Y_next)
        step = 0
        results(0, :) = [t, Y, floor%Jn]        
        if(report) then
            call update(step, n_steps, t, Y, floor)
        end if
        do while( step<n_steps )
            h_next = min(h, t_end-t)
            Y_next = rb%integrate(t, Y, h_next)
            Y = floor%contact(rb, Y_next)
            t = t + h_next
            step = step + 1
            if(report) then
                call update(step, n_steps, t, Y, floor)
            end if
            results(step, :) = [t, Y, floor%Jn]
        end do        
    end function
    
    subroutine mbd_solver(n_steps)
    integer, intent(in) :: n_steps
    
    ! Cylinder: Length=90 [mm], Radius=30 [mm]
    real(wp), parameter :: L = 90*mm
    real(wp), parameter :: D = 60*mm
    
    type(vector3) :: pos, vee, omg
    type(quaternion) :: ori
    integer :: i, n
    real(wp), allocatable :: Y(:,:)
    real(wp) :: Yi(state_size), t_end, h, omg_0
    
    integer(8) :: tic, toc, cpurate
    real(real32) :: cpu, ops
    character(1) :: key
    integer, parameter :: ESC = 27

    class(rigid_body), allocatable :: body
        
    call RANDOM_SEED()
        
    body = rigid_body( shapes(cylinder, [D/2, L]), steel_density )
    
    pos = (L/2)*k_                              ! initial position [m]
    ori = q_eye                                 ! initial orientation (quaternion)
    omg = 5.0_wp*j_                             ! initial rotational velocity [rad/s]
    vee = (-1.0_wp)*k_                          ! initial velocity [m/s]
    
    call body%set_pose(pos, ori)
    call body%set_motion(vee, omg)
        
    omg_0 = norm2(body%initial_omega)
    if( omg_0 > tiny ) then
        h = (0.5_wp*deg)/omg_0
    else
        h = 1.0_wp/n_steps
    end if
    t_end = n_steps * h    
            
    print (fmth), &
        "time", &
        "pos-x", "pos-y","pos-z", &
        "ori-x", "ori-y", "ori-z", "ori-w", &
        "mom-x", "mom-y", "mom-z", &
        "ang-x", "ang-y", "ang-z", "J"

    print (fmth), &
        "[s]", &
        "[m]", "","", &
        "[]", "", "", "", &
        "[kg m/s]", "", "", &
        "[g m^2/s]", "", "", "N*s"
    print (fmth), &
        "---", &
        "---", "---","---", &
        "---", "---", "---", "---", &
        "---", "---", "---", &
        "---", "---", "---", "---"
        
    call SYSTEM_CLOCK(tic,cpurate)    
    Y = body%simulate(t_end, n_steps, sim_update)    
    call SYSTEM_CLOCK(toc,cpurate)
            
    cpu = real(toc-tic)/cpurate
    ops = n_steps / cpu
    
    print *, ""
    print '(1x,"steps=",g0,1x,"time=",g0,1x,"kops=",g0)', n_steps, cpu, ops/1000
    
    
    end subroutine
    
    subroutine sim_update(i,n,t,Y,ctx)
    integer, intent(in) :: i,n
    real(wp), intent(in) :: t, Y(:)
    type(contact_plane), intent(in) :: ctx
        if( i==n .or. mod(i, n/report_count) ==0 ) then
            write (*, fmtv) t, Y(1:3), Y(4:7), Y(8:10), 1000*Y(11:13), ctx%Jn
        end if    
    end subroutine
    
    end module