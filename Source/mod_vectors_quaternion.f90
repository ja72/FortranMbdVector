    module mod_vectors_quaternion
    use mod_vectors_matrix3
    implicit none
    
    type :: quaternion
        type(vector3) :: vector
        real(wp) :: scalar
    contains
        procedure :: to_array => q_to_array
    end type
    
    interface assignment (=)
    module procedure :: asgn_q_array, asgn_array_q
    end interface

    interface quaternion
    module procedure :: q_from_array, q_from_vector
    end interface
    
    interface real
        module procedure :: q_to_array
    end interface    

    interface operator (+)
    module procedure :: q_add
    end interface
    interface operator (-)
    module procedure :: q_neg
    module procedure :: q_sub
    end interface
    interface operator (*)
    module procedure :: q_scale
    module procedure :: q_scale2
    end interface
    interface operator (/)
    module procedure :: q_div
    end interface
    interface dot
    module procedure :: q_dot
    end interface
    interface cross
    module procedure :: q_cross
    end interface
    interface operator (.x.)
    module procedure :: q_cross
    end interface
    
    interface operator ( .o. )
        module procedure q_product
    end interface

    interface rot
        module procedure q_axis_angle, q_vector_angle
        module procedure q_rotation_matrix
        module procedure q_rotation_zyx
    end interface
    
    interface norm2
        module procedure :: q_magnitude
    end interface
    
    interface unit
        module procedure :: q_normalize
    end interface
    
    type(quaternion), parameter :: q_eye = quaternion(o_, 1.0_wp)
        
    contains
        
    pure subroutine asgn_q_array(q,a)
    type(quaternion), intent(out) :: q
    real(wp), intent(in) :: a(4)
    q = q_from_array(a)
    end subroutine

    pure subroutine asgn_array_q(a,q)
    real(wp), intent(out) :: a(4)
    type(quaternion), intent(in) :: q
    a = q_to_array(q)
    end subroutine
    
    pure function q_from_array(a) result(q)
    real(wp), intent(in) :: a(4)
    type(quaternion) :: q
    q = quaternion(vector3(a(1),a(2),a(3)),a(4))
    end function
    
    pure function q_to_array(q) result(a)
    class(quaternion), intent(in) :: q
    real(wp) :: a(4)
        a = [q%vector%x, q%vector%y, q%vector%z, q%scalar]
    end function
    
    pure function q_from_vector(v) result(q)
    type(vector3), intent(in) :: v
    type(quaternion) :: q
        q = quaternion(v, 0.0_wp)
    end function
    
    pure function q_axis_angle(axis,angle) result(q)
    integer, intent(in) :: axis    
    real(wp), intent(in) :: angle
    type(quaternion) :: q
        q = q_vector_angle(vector3(axis),angle)
    end function
    
    pure function q_vector_angle(axis,angle) result(q)
    type(vector3), intent(in) :: axis    
    real(wp), intent(in) :: angle
    type(quaternion) :: q
    real(wp) :: s, c, m
        if(angle /= 0.0_wp) then
            m = norm2(axis)
            s = sin(angle/2)
            c = cos(angle/2)
            q = quaternion(s*axis/m,c)
        else
            q = q_eye
        end if
    end function
    
    pure function q_rotation_zyx(yaw,pitch,roll) result(q)
    ! Define quaternion rotation from axis angle
    real(wp), intent(in) :: yaw,pitch,roll
    type(quaternion) :: q
    q = rot(z_axis, yaw) .o. rot(y_axis, pitch) .o. rot(x_axis, roll)
    end function
    
    
    pure function q_magnitude(a) result(s)
        type(quaternion), intent(in) :: a
        real(wp) :: s
        s = sqrt( a%vector%x**2 + a%vector%y**2 + a%vector%z**2 + a%scalar**2 )
    end function
    
    pure function q_normalize(a) result(q)
    type(quaternion), intent(in) :: a
    type(quaternion) :: q
    real(wp) :: m
        m = norm2(a)
        if( m >= tiny ) then
            q = quaternion(a%vector/m, a%scalar/m)
        else
            q = a
        end if
    end function
    
    
    
    !-- QUATERNION ALGEBRA
    
    pure function q_add(a,b) result(c)
    type(quaternion), intent(in) :: a,b
    type(quaternion) :: c
    c%vector = a%vector + b%vector
    c%scalar = a%scalar + b%scalar
    end function

    pure function q_neg(a) result(c)
    type(quaternion), intent(in) :: a
    type(quaternion) :: c
    c%vector = -a%vector
    c%scalar = -a%scalar
    end function

    pure function q_sub(a,b) result(c)
    type(quaternion), intent(in) :: a,b
    type(quaternion) :: c
    c%vector = a%vector - b%vector
    c%scalar = a%scalar - b%scalar
    end function

    pure function q_scale(a,b) result(c)
    real(wp), intent(in) :: a
    type(quaternion), intent(in) :: b
    type(quaternion) :: c
    c%vector = a*b%vector
    c%scalar = a*b%scalar
    end function

    pure function q_scale2(b,a) result(c)
    type(quaternion), intent(in) :: b
    real(wp), intent(in) :: a
    type(quaternion) :: c
    c%vector = a*b%vector
    c%scalar = a*b%scalar
    end function

    pure function q_div(a,b) result(c)
    type(quaternion), intent(in) :: a
    real(wp), intent(in) :: b
    type(quaternion) :: c
    c%vector = a%vector/b
    c%scalar = a%scalar/b
    end function    
        
    pure function q_dot(a,b) result(s)
        type(quaternion), intent(in) :: a,b
        real(wp) :: s
        s = dot(a%vector, b%vector) + a%scalar*b%scalar
    end function
    
    pure function q_cross(a,b) result(q)
        type(quaternion), intent(in) :: a,b
        type(quaternion) :: q
        q%vector = cross(a%vector, b%vector)
        q%scalar = 0.0_wp
    end function
    
    pure function q_product(a,b) result(q)
        type(quaternion), intent(in) :: a,b
        type(quaternion) :: q
        real(wp) :: ax,ay,az,aw
        real(wp) :: bx,by,bz,bw
        !q%vector = a%scalar*b%vector + b%scalar*a%vector - cross(a%vector,b%vector)
        !q%scalar = a%scalar*b%scalar - dot(a%vector,b%vector)
        ax = a%vector%x
        ay = a%vector%y
        az = a%vector%z
        aw = a%scalar
        bx = b%vector%x
        by = b%vector%y
        bz = b%vector%z
        bw = b%scalar
        !DIR$ FMA
        q%vector = vector3( &
            aw*bx+ax*bw+ay*bz-az*by, &
            aw*by-ax*bz+ay*bw+az*bx, &
            aw*bz+ax*by-ay*bx+az*bw)
        q%scalar = aw*bw-ax*bx-ay*by-az*bz
            
    end function
    
    pure function q_rotation_matrix(q) result(R)
    ! Extract rotation matrix R from a quaternion
    ! Assumes orientation is a unit quaternion
    type(quaternion), intent(in) :: q
    type(matrix3) :: R !, X, XX    
    real(wp) :: x,y,z,w,xx,yy,zz
    
    !tex: Rotation matrix from quaternion $\boldsymbol{q}=[\vec{q_v},\,q_s]$
    !$$ \mathbf{R} = \mathbf{1} + 2\, q_s\, [\vec{q_v}\times] + 2\, [\vec{q_v}\times][\vec{q_v}\times]  $$
    
    
    !X = cross(q%vector)
    !XX = mmoi(q%vector, .true.)
    !R = eye_ + 2.0_wp*( q%scalar*X + XX )
    
    x = q%vector%x
    y = q%vector%y
    z = q%vector%z
    w = q%scalar
    xx = x*x
    yy = y*y
    zz = z*z
    
    !DIR$ FMA
    R%a11 = 1 - 2*(yy+zz)
    R%a22 = 1 - 2*(xx+zz)
    R%a33 = 1 - 2*(xx+yy)
    R%a12 = -2*(w*z-x*y)
    R%a21 = -2*(-w*z-x*y)
    R%a13 = -2*(-x*z-w*y)
    R%a31 = -2*(-x*z+w*y)
    R%a23 = -2*( w*x-y*z)
    R%a32 = -2*(-w*x-y*z)
    
    end function
    
    end module