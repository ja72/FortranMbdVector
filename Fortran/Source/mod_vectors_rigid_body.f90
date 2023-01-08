    module mod_vectors_rigid_body
    use mod_vectors_quaternion
    use mod_vectors_shapes
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
        rb%initial_position = vector3(0.0_wp,0.0_wp,0.0_wp)
        rb%initial_orientation = quaternion(vector3(0.0_wp,0.0_wp,0.0_wp), 0.0_wp)
        rb%initial_velocity = vector3(0.0_wp,0.0_wp,0.0_wp)
        rb%initial_omega = vector3(0.0_wp,0.0_wp,0.0_wp)
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
        
        
    
    end module