    module mod_vectors_contact
    use mod_vectors_rigid_body
    implicit none
    
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
    
    
    
    contains
    
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
        ctx%slip   = vector3(0.0_wp, 0.0_wp, 0.0_wp)
        ctx%delta  = 0.0_wp
        ctx%v_imp  = 0.0_wp
        ctx%v_slip = 0.0_wp
        ctx%Jn     = 0.0_wp
        ctx%Js     = 0.0_wp
        ctx%active = .false.
    end function
    
    pure subroutine ctx_reset(ctx)
    class(contact_plane), intent(inout) :: ctx
        ctx%slip   = vector3(0.0_wp, 0.0_wp, 0.0_wp)
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
    type(matrix3) :: R, Rt, I_inv, cx, M_inv
    real(wp) :: m_eff, m_slip
            
        ! Body motions
        pos_b = Y(1:3)
        ori = Y(4:7) 
        p = Y(8:10)
        L_b = Y(11:13)
        R = rot(ori)
        Rt = rot(ori, inverse = .true.)
        c = R * rb%cg
        I_inv = R*diag(1/rb%mmoi)*Rt
        omg = I_inv*(L_b - cross(c, p))
        vee_b = p/rb%mass - cross(omg, c)
        
        ! contact point on surface of shape
        r_A = pos_b + c + R*rb%shape%nearest_point(Rt*ctx%normal)
        ctx%delta = dot(ctx%normal, r_A - ctx%origin)
        ! clip contact point on plane
        r_A = r_A - min(ctx%delta, 0.0_wp)*ctx%normal
                
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
                ctx%slip = vector3(0.0_wp, 0.0_wp, 0.0_wp)
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
            ctx%slip   = vector3(0.0_wp, 0.0_wp, 0.0_wp)
            ctx%v_slip = 0.0_wp
            ctx%Jn     = 0.0_wp
            ctx%Js     = 0.0_wp
            ! Return body state unmodified
            Y_next = Y
        end if
    end function        
    

    end module