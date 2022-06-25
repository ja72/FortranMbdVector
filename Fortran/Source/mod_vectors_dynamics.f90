    module mod_vectors_dynamics
    use mod_vectors_quaternion
    use mod_vectors_geometry
    implicit none
        
    ! Simulation parameters
    integer, parameter :: state_size = 13
        
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
        procedure :: est_max_time_step => rb_get_max_time_step
    end type
    
    interface rigid_body
        module procedure :: rb_from_shape, rb_new
    end interface
    
    ! Contact Surface Definition
    type :: contact_plane
        real(wp) :: epsilon, mu
        type(vector3) :: origin
        type(vector3) :: normal
        type(vector3) :: slip
        real(wp) :: delta, v_imp, Jn
        real(wp) :: v_slip, Js
        logical :: active
    contains
        procedure :: calculate => ctx_calculate
        procedure :: reset => ctx_reset
    end type
        
    interface contact_plane
        module procedure :: ctx_new_plane
    end interface
    
    ! Simulation contains 1 body and 1 contact
    type :: simulation
        type(rigid_body) :: body
        type(contact_plane) :: contact
        type(vector3) :: gravity
        integer :: step
        real(wp) :: time
        real(wp), dimension(state_size) :: current
        procedure(step_handler), pointer :: step_taken
    contains
        procedure :: reset => sim_reset
        procedure :: run => sim_run
        procedure :: integrate => sim_integrate
    end type
    
    interface simulation
        module procedure :: new_simulation
    end interface
    
    interface 
        subroutine step_handler(sim,n)
        import
        integer, intent(in) :: n
        class(simulation), intent(in) :: sim
        end subroutine
    end interface
            
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
    
    function rb_from_shape(s, density, mass) result(rb)
    class(shapes), intent(in) :: s
    real(wp), intent(in), optional :: density
    real(wp), intent(in), optional :: mass
    type(rigid_body) :: rb
    real(wp) :: m

    if( present(mass) ) then
        m = mass
    elseif ( present(density) ) then
        m = density * s%volume
    else
        m = 1000.0_wp * s%volume
    end if
    rb = rigid_body(s, mass, s%center, mass * s%vmmoi)

    end function
    
    pure function rb_get_max_time_step(rb, n_steps, end_time, angularResolutionDegrees) result(h)
    ! Find maximum time step based on recommended angular resolution for rotating
    ! bodies. Defaults to 1 deg of rotation per frame.
    class(rigid_body), intent(in) :: rb
    integer, intent(in) :: n_steps
    real(wp), intent(in), optional :: end_time, angularResolutionDegrees
    real(wp) :: h, omg    
        if(present(end_time)) then
            h = end_time / n_steps
        else
            h = 1.0_wp / n_steps
        end if
        omg = norm2(rb%initial_velocity)
        if( omg > tiny) then
            if(present(angularResolutionDegrees)) then
                h = min(h, angularResolutionDegrees * deg / omg)
            else
                h = min(h, 1.0_wp * deg / omg)
            end if
        end if    
    end function
    
    pure subroutine rb_set_cg_pose(rb, pos, ori)
    class(rigid_body), intent(inout) :: rb
    type(vector3), intent(in) :: pos
    type(quaternion), intent(in) :: ori
    type(vector3) :: c
        c = rot(ori, rb%cg)
        rb%initial_position = pos - c
        rb%initial_orientation = ori
    end subroutine
    
    pure subroutine rb_set_cg_motion(rb, vee, omg)
    class(rigid_body), intent(inout) :: rb
    type(vector3), intent(in) :: vee, omg
    type(vector3) :: c
        c = rot(rb%initial_orientation, rb%cg)
        rb%initial_velocity = vee + cross(c,omg)
        rb%initial_omega = omg
    end subroutine
        
    pure function rb_get_state(rb, pos_b, ori, vee_b, omg) result(Y)
    ! Calculate body state from motion vectors
    class(rigid_body), intent(in) :: rb
    type(vector3), intent(in) :: pos_b, vee_b, omg
    type(quaternion), intent(in) :: ori
    real(wp):: Y(state_size)
    type(vector3) :: p, L_b, c
    type(matrix3) :: R, I_c
        !tex: Get momentum at reference point from motion
            ! $$\begin{aligned}\vec{p} & =m\left(\vec{v}_{b}+\vec{\omega}\times\vec{c}\right)\\
            !\vec{L}_{b} & ={\rm I}_{c}\vec{\omega}+\vec{c}\times\vec{p}
            !\end{aligned}$$    
        R = rot(ori)
        c = R * rb%cg
        I_c = R*diag(rb%mmoi)*R%transpose()
        p = rb%mass*(vee_b - cross(c, omg))
        L_b = I_c*omg + cross(c, p)        
        Y = [ &
            pos_b%x, pos_b%y, pos_b%z, &
            ori%vector%x, ori%vector%y, ori%vector%z, ori%scalar, &
            p%x, p%y, p%z, &
            L_b%x, L_b%y, L_b%z]
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
    
    pure subroutine rb_get_motion(rb, Y, vee_b, omg)
    ! Calculate motion vectors from body state
    class(rigid_body), intent(in) :: rb
    real(wp), intent(in) :: Y(state_size)
    type(vector3), intent(out) :: vee_b
    type(vector3), intent(out) :: omg
    type(vector3) :: c, p, L_b
    type(quaternion) :: ori
    type(matrix3) :: R, I_inv
        !tex: Get motion from momentum
        ! $$\begin{aligned}\vec{v}_{b} & =\tfrac{1}{m}\vec{p}-\vec{\omega}\times\vec{c}\\
        !\vec{\omega} & ={\rm I}_{c}^{-1}\left(\vec{L}_{b}-\vec{c}\times\vec{p}\right)
        !\end{aligned}$$
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
    
    pure function rb_get_rate(rb, t, Y, gee) result(Yp)
    ! Calculate the body state time derivative Yp=f(t,Y)
    ! NOTE: Needs input from simulation conditions for
    ! gravity field (gee) for weight calculation.
    class(rigid_body), intent(in) :: rb
    real(wp), intent(in) :: t, Y(state_size)
    type(vector3), intent(in) :: gee
    real(wp) :: Yp(state_size)
    type(vector3) :: F, tau_b, omg, vee_b
    type(quaternion) :: ori, qp
    type(vector3) :: c, p, L_b, dLdt
    type(matrix3) :: R, I_inv
    
        call y_normalize(Y)
    
        ori = Y(4:7) 
        p = Y(8:10)
        L_b = Y(11:13)
        R = rot(ori)
        c = R * rb%cg
        I_inv = R*diag(1/rb%mmoi)*R%transpose()
        omg = I_inv*(L_b - cross(c, p))
        vee_b = p/rb%mass - cross(omg, c)
        qp = 0.5_wp*(quaternion(omg) .o. ori)
        F = rb%mass * gee
        tau_b = cross(c, F)
        dLdt = tau_b - cross(vee_b, p)
        Yp = [ &
            vee_b%x, vee_b%y, vee_b%z, &
            qp%vector%x, qp%vector%y, qp%vector%z, qp%scalar, &
            F%x, F%y, F%z, &
            dLdt%x, dLdt%y, dLdt%z]
    
    end function
        
    pure function ctx_new_plane(origin,normal,epsilon,mu) result(ctx)
    type(contact_plane) :: ctx
    real(wp),intent(in), optional :: epsilon, mu
    type(vector3),intent(in) :: origin
    type(vector3),intent(in) :: normal    
    
        if(present(epsilon)) then
            ctx%epsilon = epsilon
        else
            ctx%epsilon = 0.0_wp
        endif
        
        if(present(mu)) then
            ctx%mu = mu
        else
            ctx%mu = 0.0_wp
        endif
        
        ctx%origin = origin
        ctx%normal = normal/norm2(normal)
        ctx%slip   = o_
        ctx%delta  = 0.0_wp
        ctx%v_imp  = 0.0_wp
        ctx%v_slip = 0.0_wp
        ctx%Jn     = 0.0_wp
        ctx%Js     = 0.0_wp
        ctx%active = .false.
    end function
    
    pure subroutine ctx_reset(ctx)
    class(contact_plane), intent(inout) :: ctx
        ctx%slip   = o_
        ctx%delta  = 0.0_wp
        ctx%v_imp  = 0.0_wp
        ctx%v_slip = 0.0_wp
        ctx%Jn     = 0.0_wp
        ctx%Js     = 0.0_wp
        ctx%active = .false.
    end subroutine
        
    function ctx_calculate(ctx, rb, Y) result(Y_next)
    class(contact_plane), intent(inout) :: ctx
    type(rigid_body), intent(in) :: rb
    real(wp), intent(in) :: Y(state_size)
    real(wp) :: Y_next(state_size)
    type(vector3) :: p, L_b, omg, vee_b, pos_b
    type(quaternion) :: ori
    type(vector3) :: c, r_A, c_A, v_A, J_imp
    type(matrix3) :: R, I_inv, cx, M_inv
    real(wp) :: m_eff, m_slip
            
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
                
        ! get contact speed
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
            
            ! Find sliding velociy properties
            ctx%slip = v_A - ctx%v_imp * ctx%normal
            ctx%v_slip = norm2(ctx%slip)
            if( abs(ctx%v_slip) < tiny ) then
                ! negligible sliding velocity
                ! NOTE: sticking isn't implemented yet
                ctx%slip = o_
                ctx%v_slip = 0.0_wp
                ctx%Js = 0.0_wp
            else
                ! non-neglible sliding velocity
                ctx%slip = -ctx%slip/ctx%v_slip
                m_slip = 1/dot(ctx%slip, M_inv*ctx%slip)
                ! find sliding impulse to stop the movement
                ctx%Js = m_slip*ctx%v_slip
                ! cap sliding impulse to traction limit
                ctx%Js = min(ctx%mu*ctx%Jn, ctx%Js)
            end if
                        
            ! Combined impulse included normal and sliding components
            J_imp = ctx%Jn*ctx%normal + ctx%Js*ctx%slip
            p = p + J_imp
            L_b = L_b + cross(r_A-pos_b, J_imp)
            ! Return new body state with upadted momenta
            Y_next = [ &
                pos_b%x, pos_b%y, pos_b%z, &
                ori%vector%x, ori%vector%y, ori%vector%z, ori%scalar, &
                p%x, p%y, p%z, &
                L_b%x, L_b%y, L_b%z]
        else
            ! Contat not active
            ctx%slip   = o_
            ctx%v_slip = 0.0_wp
            ctx%Jn     = 0.0_wp
            ctx%Js     = 0.0_wp
            ! Return body state unmodified
            Y_next = Y
        end if
    end function        
        
    function new_simulation(rb,ctx,gee,update) result(sim)
    type(simulation) :: sim
    type(rigid_body), intent(in) :: rb
    type(contact_plane), intent(in) :: ctx
    type(vector3), intent(in) :: gee
    procedure(step_handler), pointer, intent(in), optional :: update
        sim%body = rb
        sim%contact = ctx
        sim%gravity = gee
        if(present(update)) then
            sim%step_taken => update
        else
            sim%step_taken => null()
        end if
        call sim%reset()
    
    end function
    
    subroutine sim_reset(sim)
    class(simulation), intent(inout) :: sim
        sim%step = 0
        sim%time = 0.0_wp
        sim%current = sim%body%initial_state()
        sim%current = sim%contact%calculate(sim%body, sim%current)
    end subroutine
    
    pure function sim_integrate(sim, h) result(Y_next)
    ! Implement RK4 integrator for rigid bodies
    class(simulation), intent(in) :: sim
    real(wp), intent(in) :: h
    real(wp), dimension(state_size) :: Y_next, Y1, Y2, Y3
    real(wp) :: k0(state_size),k1(state_size),k2(state_size),k3(state_size)
    
        !DEC$ FMA
        k0 = rb_get_rate(sim%body, sim%time    , sim%current, sim%gravity)
        Y1 = sim%current + (h/2)*k0
        k1 = rb_get_rate(sim%body, sim%time+h/2, Y1, sim%gravity)
        Y2 = sim%current + (h/2)*k1
        k2 = rb_get_rate(sim%body, sim%time+h/2, Y2, sim%gravity)
        Y3 = sim%current + (h)*k2
        k3 = rb_get_rate(sim%body, sim%time+h  , Y3, sim%gravity)
        
        Y_next = sim%current + (h/6)*k0+(h/3)*k1+(h/3)*k2+(h/6)*k3
        call y_normalize(Y_next)
    end function
    
    subroutine sim_run(sim, t_end, n_steps)
    class(simulation), intent(inout) :: sim
    integer, intent(in) :: n_steps
    real(wp), intent(in) :: t_end
    real(wp) :: h, h_next
    real(wp), dimension(state_size) :: next
    integer :: n
        ! loop integrator and contact until end time is reached
        if( associated(sim%step_taken)) then
            call sim%step_taken(n_steps)
        end if
        h = (t_end - sim%time) / n_steps
        n = n_steps + sim%step
        do while(sim%step < n)
            h_next = min(h, t_end - sim%time)
            next = sim_integrate(sim, h_next)
            sim%current = ctx_calculate(sim%contact, sim%body, next)
            sim%time = sim%time + h_next
            sim%step = sim%step + 1
            if( associated(sim%step_taken)) then
                call sim%step_taken(n_steps)
            end if
        end do
    
    end subroutine
    
    end module