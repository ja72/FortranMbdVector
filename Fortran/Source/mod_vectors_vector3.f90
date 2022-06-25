    module mod_vectors_vector3
    use mod_constants
    implicit none
    
    ! Common directions
    enum, bind(c)
        enumerator :: x_axis = 1
        enumerator :: y_axis = 2
        enumerator :: z_axis = 3
    end enum
    
    type :: vector3
        real(wp) :: x,y,z
    contains
        procedure :: to_array => vec3_to_array
    end type

    interface assignment (=)
    module procedure :: asgn_vec3_array, asgn_array_vec3
    end interface

    interface vector3
    module procedure :: vec3_from_array, vec3_from_axis
    end interface
    
    interface norm2
        module procedure :: vec3_magnitude
    end interface
    
    interface real
        module procedure :: vec3_to_array
    end interface

    interface operator (+)
    module procedure :: vec3_add
    end interface
    interface operator (-)
    module procedure :: vec3_neg
    module procedure :: vec3_sub
    end interface
    interface operator (*)
    module procedure :: vec3_scale
    module procedure :: vec3_scale2
    end interface
    interface operator (/)
    module procedure :: vec3_div
    end interface
    interface dot
    module procedure :: vec3_dot
    end interface
    interface cross
    module procedure :: vec3_cross
    end interface
    interface operator (.x.)
    module procedure :: vec3_cross
    end interface

    type(vector3), parameter :: o_ = vector3(0.0_wp,0.0_wp,0.0_wp)
    type(vector3), parameter :: i_ = vector3(1.0_wp,0.0_wp,0.0_wp)
    type(vector3), parameter :: j_ = vector3(0.0_wp,1.0_wp,0.0_wp)
    type(vector3), parameter :: k_ = vector3(0.0_wp,0.0_wp,1.0_wp)

    contains

    pure subroutine asgn_vec3_array(v,a)
    type(vector3), intent(out) :: v
    real(wp), intent(in) :: a(3)
    v = vec3_from_array(a)
    end subroutine

    pure subroutine asgn_array_vec3(a,v)
    real(wp), intent(out) :: a(3)
    type(vector3), intent(in) :: v
    a = vec3_to_array(v)
    end subroutine
    
    pure function vec3_from_axis(axis) result(v)
    integer, intent(in) :: axis    
    type(vector3) :: v
        v%x = 0.0_wp
        v%y = 0.0_wp
        v%z = 0.0_wp
        select case(axis)
        case (x_axis)
            v%x = 1.0_wp
        case (y_axis)
            v%y = 1.0_wp
        case (z_axis)
            v%z = 1.0_wp
        end select
    end function
    
    pure function vec3_from_array(a) result(v)
    real(wp), intent(in) :: a(3)
    type(vector3) :: v
    v%x = a(1)
    v%y = a(2)
    v%z = a(3)
    end function

    pure function vec3_to_array(v) result(a)
    real(wp) :: a(3)
    class(vector3), intent(in) :: v
    a = [v%x, v%y, v%z]
    end function

    pure function vec3_add(a,b) result(c)
    type(vector3), intent(in) :: a,b
    type(vector3) :: c
    c%x = a%x + b%x
    c%y = a%y + b%y
    c%z = a%z + b%z
    end function

    pure function vec3_neg(a) result(c)
    type(vector3), intent(in) :: a
    type(vector3) :: c
    c%x = -a%x
    c%y = -a%y
    c%z = -a%z
    end function

    pure function vec3_sub(a,b) result(c)
    type(vector3), intent(in) :: a,b
    type(vector3) :: c
    c%x = a%x - b%x
    c%y = a%y - b%y
    c%z = a%z - b%z
    end function

    pure function vec3_scale(a,b) result(c)
    real(wp), intent(in) :: a
    type(vector3), intent(in) :: b
    type(vector3) :: c
    c%x = a*b%x
    c%y = a*b%y
    c%z = a*b%z
    end function

    pure function vec3_scale2(b,a) result(c)
    type(vector3), intent(in) :: b
    real(wp), intent(in) :: a
    type(vector3) :: c
    c%x = a*b%x
    c%y = a*b%y
    c%z = a*b%z
    end function

    pure function vec3_div(a,b) result(c)
    type(vector3), intent(in) :: a
    real(wp), intent(in) :: b
    type(vector3) :: c
    c%x = a%x/b
    c%y = a%y/b
    c%z = a%z/b
    end function
    
    pure function vec3_magnitude(v) result(s)
    type(vector3), intent(in) :: v
    real(wp) :: s
        s = sqrt(v%x**2 + v%y**2 + v%z**2 )
    end function
    
    pure function vec3_dot(a,b) result(s)
    type(vector3), intent(in) :: a,b
    real(wp) :: s
        s = a%x*b%x + a%y*b%y + a%z*b%z
    end function
    
    pure function vec3_cross(a,b) result(v)
    type(vector3), intent(in) :: a,b
    type(vector3) :: v
        v%x = a%y*b%z - a%z*b%y
        v%y = a%z*b%x - a%x*b%z
        v%z = a%x*b%y - a%y*b%x
    end function
    
    end module