    module mod_vectors_matrix3
    use mod_vectors_vector3
    implicit none

    type :: matrix3
        real(wp) :: a11,a21,a31, a12,a22,a32, a13,a23,a33
    contains
    procedure :: to_array => mat3_to_array
    procedure :: to_matrix => mat3_to_matrix
    procedure :: symm => mat3_to_smat3
    procedure :: column => mat3_get_column
    procedure :: row => mat3_get_row
    procedure :: transpose => mat3_transpose
    end type

    type :: symmatrix3
        real(wp) :: a11, a21, a31, a22, a32, a33
    contains
    procedure :: to_array => smat3_to_array
    procedure :: to_matrix => smat3_to_matrix
    procedure :: full => smat3_to_mat3
    procedure :: column => smat3_get_column
    procedure :: row => smat3_get_row
    procedure :: transpose => smat3_transpose
    end type
    
    interface real
        module procedure :: mat3_to_array, smat3_to_array
    end interface
    

    interface assignment (=)
    module procedure :: asgn_mat3_array, asgn_array_mat3
    module procedure :: asgn_smat3_array, asgn_array_smat3
    module procedure :: asgn_mat3_matrix, asgn_matrix_mat3
    module procedure :: asgn_smat3_matrix, asgn_matrix_smat3
    module procedure :: asgn_mat3_smat3
    end interface

    interface matrix3
    module procedure :: mat3_from_array, mat3_from_matrix
    module procedure :: smat3_to_mat3
    end interface
    interface symmatrix3
    module procedure :: smat3_from_array, smat3_from_matrix
    module procedure :: mat3_to_smat3
    end interface

    interface operator (+)
    module procedure :: mat3_add, smat3_add
    end interface
    interface operator (-)
    module procedure :: mat3_neg, smat3_neg
    module procedure :: mat3_sub, smat3_sub
    end interface
    interface operator (*)
    module procedure :: mat3_scale, smat3_scale
    module procedure :: mat3_scale2, smat3_scale2
    module procedure :: mat3_product_vec3, vec3_product_mat3
    module procedure :: smat3_product_vec3, vec3_product_smat3
    module procedure :: mat3_product_mat3, mat3_product_smat3
    module procedure :: smat3_product_mat3, smat3_product_smat3
    end interface
    interface operator (/)
    module procedure :: mat3_div, smat3_div
    end interface
    interface diag
    module procedure :: smat3_from_diag
    end interface
    interface cross
    module procedure :: vec3_cross_matrix
    end interface
    interface outer
    module procedure :: vec3_outer
    end interface
    interface mmoi
    module procedure :: vec3_mmoi
    end interface

    type(matrix3), parameter :: zeros_ = matrix3(0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp)
    type(matrix3), parameter :: eye_ = matrix3(1.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp)

    contains

    pure subroutine asgn_mat3_array(m,a)
    type(matrix3), intent(out) :: m
    real(wp), intent(in) :: a(9)
    m = mat3_from_array(a)
    end subroutine

    pure subroutine asgn_mat3_matrix(m,a)
    type(matrix3), intent(out) :: m
    real(wp), intent(in) :: a(3,3)
    m = mat3_from_matrix(a)
    end subroutine

    pure subroutine asgn_array_mat3(a,m)
    real(wp), intent(out) :: a(9)
    type(matrix3), intent(in) :: m
    a = mat3_to_array(m)
    end subroutine
    
    pure subroutine asgn_matrix_mat3(a,m)
    real(wp), intent(out) :: a(3,3)
    type(matrix3), intent(in) :: m
    a = mat3_to_matrix(m)
    end subroutine

    pure subroutine asgn_smat3_array(m,a)
    type(symmatrix3), intent(out) :: m
    real(wp), intent(in) :: a(6)
    m = smat3_from_array(a)
    end subroutine
    
    pure subroutine asgn_smat3_matrix(m,a)
    type(symmatrix3), intent(out) :: m
    real(wp), intent(in) :: a(3,3)
    m = smat3_from_matrix(a)
    end subroutine

    pure subroutine asgn_array_smat3(a,m)
    real(wp), intent(out) :: a(6)
    type(symmatrix3), intent(in) :: m
    a = smat3_to_array(m)
    end subroutine

    pure subroutine asgn_mat3_smat3(m,s)
    type(matrix3), intent(out) :: m
    type(symmatrix3), intent(in) :: s
    m = smat3_to_matrix(s)
    end subroutine

    pure subroutine asgn_matrix_smat3(a,m)
    real(wp), intent(out) :: a(3,3)
    type(symmatrix3), intent(in) :: m
    a = smat3_to_matrix(m)
    end subroutine

    pure function mat3_from_matrix(a) result(m)
    real(wp), intent(in) :: a(3,3)
    type(matrix3) :: m
    m = mat3_from_array( reshape(a, [9]) )
    end function

    pure function mat3_from_array(a) result(m)
    real(wp), intent(in) :: a(9)
    type(matrix3) :: m
    m = matrix3(a(1), a(2), a(3), &
        a(4), a(5), a(6), &
        a(7), a(8), a(9))
    end function

    pure function smat3_from_matrix(a) result(m)
    real(wp), intent(in) :: a(3,3)
    type(symmatrix3) :: m
    m = symmatrix3(a(1,1), (a(2,1)+a(1,2))/2, &
        (a(3,1)+a(3,1))/2, a(2,2), &
        (a(2,3)+a(3,2))/2, a(3,3) )
    end function

    pure function smat3_from_array(a) result(m)
    real(wp), intent(in) :: a(6)
    type(symmatrix3) :: m
    m = symmatrix3(a(1),a(2),a(3),a(4),a(5),a(6) )
    end function

    pure function mat3_to_array(m) result(a)
    real(wp) :: a(9)
    class(matrix3), intent(in) :: m
    a = [m%a11, m%a21, m%a31, &
        m%a12, m%a22, m%a32, &
        m%a13, m%a23, m%a33]
    end function

    pure function mat3_to_matrix(m) result(a)
    real(wp) :: a(3,3)
    class(matrix3), intent(in) :: m
    a = reshape( [m%a11, m%a21, m%a31, m%a12, m%a22, m%a32, m%a13, m%a23, m%a33], [3,3] )
    end function

    pure function mat3_to_smat3(m) result(a)
    class(matrix3), intent(in) :: m
    type(symmatrix3) :: a
    a = symmatrix3( [ m%a11, (m%a12+m%a21)/2, (m%a13+m%a31)/2, &
        m%a22, (m%a23+m%a32)/2, m%a33] )
    end function
    
    pure function smat3_to_mat3(m) result(a)
    class(symmatrix3), intent(in) :: m
    type(matrix3) :: a
    a%a11 = m%a11
    a%a22 = m%a22
    a%a33 = m%a33
    a%a12 = m%a21
    a%a13 = m%a31
    a%a23 = m%a32
    a%a21 = m%a21
    a%a31 = m%a31
    a%a32 = m%a32
    end function

    pure function smat3_to_array(m) result(a)
    real(wp) :: a(6)
    class(symmatrix3), intent(in) :: m
    a = [m%a11, m%a21, m%a31, m%a22, m%a32, m%a33]
    end function

    pure function smat3_to_matrix(m) result(a)
    real(wp) :: a(3,3)
    class(symmatrix3), intent(in) :: m
    a = reshape( [m%a11, m%a21, m%a31, m%a21, m%a22, m%a32, m%a31, m%a32, m%a33], [3,3] )
    end function
    
    pure function smat3_from_diag(d) result(m)
    real(wp), intent(in) :: d(3)
    type(symmatrix3) :: m
        m%a11 = d(1)
        m%a22 = d(2)
        m%a33 = d(3)
        m%a21 = 0.0_wp
        m%a31 = 0.0_wp
        m%a32 = 0.0_wp
    end function
    
    
    !-- MATRIX3 ALGEBRA
    pure function mat3_add(a,b) result(c)
    type(matrix3), intent(in) :: a,b
    type(matrix3) :: c
    !DIR$ FMA
    c%a11 = a%a11 + b%a11
    c%a12 = a%a12 + b%a12
    c%a13 = a%a13 + b%a13
    c%a21 = a%a21 + b%a21
    c%a22 = a%a22 + b%a22
    c%a23 = a%a23 + b%a23
    c%a31 = a%a31 + b%a31
    c%a32 = a%a32 + b%a32
    c%a33 = a%a33 + b%a33
    end function
    pure function mat3_neg(a) result(c)
    type(matrix3), intent(in) :: a
    type(matrix3) :: c
    !DIR$ FMA
    c%a11 = -a%a11
    c%a12 = -a%a12
    c%a13 = -a%a13
    c%a21 = -a%a21
    c%a22 = -a%a22
    c%a23 = -a%a23
    c%a31 = -a%a31
    c%a32 = -a%a32
    c%a33 = -a%a33
    end function
    pure function mat3_sub(a,b) result(c)
    type(matrix3), intent(in) :: a,b
    type(matrix3) :: c
    !DIR$ FMA
    c%a11 = a%a11 - b%a11
    c%a12 = a%a12 - b%a12
    c%a13 = a%a13 - b%a13
    c%a21 = a%a21 - b%a21
    c%a22 = a%a22 - b%a22
    c%a23 = a%a23 - b%a23
    c%a31 = a%a31 - b%a31
    c%a32 = a%a32 - b%a32
    c%a33 = a%a33 - b%a33
    end function
    
    pure function mat3_scale(a,b) result(c)
    real(wp), intent(in) :: a
    type(matrix3), intent(in) :: b
    type(matrix3) :: c
    !DIR$ FMA
    c%a11 = a*b%a11
    c%a12 = a*b%a12
    c%a13 = a*b%a13
    c%a21 = a*b%a21
    c%a22 = a*b%a22
    c%a23 = a*b%a23
    c%a31 = a*b%a31
    c%a32 = a*b%a32
    c%a33 = a*b%a33
    end function

    pure function mat3_scale2(b,a) result(c)
    type(matrix3), intent(in) :: b
    real(wp), intent(in) :: a
    type(matrix3) :: c
    !DIR$ FMA
    c%a11 = a*b%a11
    c%a12 = a*b%a12
    c%a13 = a*b%a13
    c%a21 = a*b%a21
    c%a22 = a*b%a22
    c%a23 = a*b%a23
    c%a31 = a*b%a31
    c%a32 = a*b%a32
    c%a33 = a*b%a33
    end function

    pure function mat3_div(a,b) result(c)
    type(matrix3), intent(in) :: a
    real(wp), intent(in) :: b
    type(matrix3) :: c
    !DIR$ FMA
    c%a11 = a%a11/b
    c%a12 = a%a12/b
    c%a13 = a%a13/b
    c%a21 = a%a21/b
    c%a22 = a%a22/b
    c%a23 = a%a23/b
    c%a31 = a%a31/b
    c%a32 = a%a32/b
    c%a33 = a%a33/b
    end function

    pure function mat3_product_vec3(a,b) result(c)
    type(matrix3), intent(in) :: a
    type(vector3), intent(in) :: b
    type(vector3) :: c
    !DIR$ FMA
    c%x = a%a11*b%x + a%a12*b%y + a%a13*b%z
    c%y = a%a21*b%x + a%a22*b%y + a%a23*b%z
    c%z = a%a31*b%x + a%a32*b%y + a%a33*b%z
    end function

    pure function vec3_product_mat3(b,a) result(c)
    type(vector3), intent(in) :: b
    type(matrix3), intent(in) :: a
    type(vector3) :: c
    !DIR$ FMA
    c%x = a%a11*b%x + a%a21*b%y + a%a31*b%z
    c%y = a%a12*b%x + a%a22*b%y + a%a32*b%z
    c%z = a%a13*b%x + a%a23*b%y + a%a33*b%z
    end function
    pure function mat3_product_mat3(a,b) result(c)
    type(matrix3), intent(in) :: a
    type(matrix3), intent(in) :: b
    type(matrix3) :: c
    !DIR$ FMA
    c%a11 = a%a11*b%a11 + a%a12*b%a21 + a%a13*b%a31
    c%a12 = a%a11*b%a12 + a%a12*b%a22 + a%a13*b%a32
    c%a13 = a%a11*b%a13 + a%a12*b%a23 + a%a13*b%a33
    c%a21 = a%a21*b%a11 + a%a22*b%a21 + a%a23*b%a31
    c%a22 = a%a21*b%a12 + a%a22*b%a22 + a%a23*b%a32
    c%a23 = a%a21*b%a13 + a%a22*b%a23 + a%a23*b%a33
    c%a31 = a%a31*b%a11 + a%a32*b%a21 + a%a33*b%a31
    c%a32 = a%a31*b%a12 + a%a32*b%a22 + a%a33*b%a32
    c%a33 = a%a31*b%a13 + a%a32*b%a23 + a%a33*b%a33    
    end function
    
    
    !-- SYMMATRIX3 ALGEBRA
    pure function smat3_add(a,b) result(c)
    type(symmatrix3), intent(in) :: a,b
    type(symmatrix3) :: c
    !DIR$ FMA
    c%a11 = a%a11 + b%a11
    c%a21 = a%a21 + b%a21
    c%a22 = a%a22 + b%a22
    c%a31 = a%a31 + b%a31
    c%a32 = a%a32 + b%a32
    c%a33 = a%a33 + b%a33
    end function
    pure function smat3_neg(a) result(c)
    type(symmatrix3), intent(in) :: a
    type(symmatrix3) :: c
    !DIR$ FMA
    c%a11 = - a%a11
    c%a21 = - a%a21
    c%a22 = - a%a22
    c%a31 = - a%a31
    c%a32 = - a%a32
    c%a33 = - a%a33
    end function
    
    pure function smat3_sub(a,b) result(c)
    type(symmatrix3), intent(in) :: a,b
    type(symmatrix3) :: c
    !DIR$ FMA
    c%a11 = a%a11 - b%a11
    c%a21 = a%a21 - b%a21
    c%a22 = a%a22 - b%a22
    c%a31 = a%a31 - b%a31
    c%a32 = a%a32 - b%a32
    c%a33 = a%a33 - b%a33
    end function
    pure function smat3_scale(a,b) result(c)
    real(wp), intent(in) :: a
    type(symmatrix3), intent(in) :: b
    type(symmatrix3) :: c
    !DIR$ FMA
    c%a11 = a*b%a11
    c%a21 = a*b%a21
    c%a22 = a*b%a22
    c%a31 = a*b%a31
    c%a32 = a*b%a32
    c%a33 = a*b%a33
    end function
    pure function smat3_scale2(b,a) result(c)
    type(symmatrix3), intent(in) :: b
    real(wp), intent(in) :: a
    type(symmatrix3) :: c
    !DIR$ FMA
    c%a11 = a*b%a11
    c%a21 = a*b%a21
    c%a22 = a*b%a22
    c%a31 = a*b%a31
    c%a32 = a*b%a32
    c%a33 = a*b%a33
    end function
    pure function smat3_div(a,b) result(c)
    type(symmatrix3), intent(in) :: a
    real(wp), intent(in) :: b
    type(symmatrix3) :: c
    !DIR$ FMA
    c%a11 = a%a11/b
    c%a21 = a%a21/b
    c%a22 = a%a22/b
    c%a31 = a%a31/b
    c%a32 = a%a32/b
    c%a33 = a%a33/b
    end function
    pure function smat3_product_vec3(a,b) result(c)
    type(symmatrix3), intent(in) :: a
    type(vector3), intent(in) :: b
    type(vector3) :: c
    !DIR$ FMA
    c%x = a%a11*b%x + a%a21*b%y + a%a31*b%z
    c%y = a%a21*b%x + a%a22*b%y + a%a32*b%z
    c%z = a%a31*b%x + a%a32*b%y + a%a33*b%z
    end function
    pure function vec3_product_smat3(b,a) result(c)
    type(vector3), intent(in) :: b
    type(symmatrix3), intent(in) :: a
    type(vector3) :: c
    !DIR$ FMA
    c%x = a%a11*b%x + a%a21*b%y + a%a31*b%z
    c%y = a%a21*b%x + a%a22*b%y + a%a32*b%z
    c%z = a%a31*b%x + a%a32*b%y + a%a33*b%z
    end function

    pure function mat3_product_smat3(a,b) result(c)
    type(matrix3), intent(in) :: a
    type(symmatrix3), intent(in) :: b
    type(matrix3) :: c
    !DIR$ FMA
    c%a11 = a%a11*b%a11 + a%a12*b%a21 + a%a13*b%a31
    c%a12 = a%a11*b%a21 + a%a12*b%a22 + a%a13*b%a32
    c%a13 = a%a11*b%a31 + a%a12*b%a32 + a%a13*b%a33
    c%a21 = a%a21*b%a11 + a%a22*b%a21 + a%a23*b%a31
    c%a22 = a%a21*b%a21 + a%a22*b%a22 + a%a23*b%a32
    c%a23 = a%a21*b%a31 + a%a22*b%a32 + a%a23*b%a33
    c%a31 = a%a31*b%a11 + a%a32*b%a21 + a%a33*b%a31
    c%a32 = a%a31*b%a21 + a%a32*b%a22 + a%a33*b%a32
    c%a33 = a%a31*b%a31 + a%a32*b%a32 + a%a33*b%a33    
    end function
    pure function smat3_product_mat3(a,b) result(c)
    type(symmatrix3), intent(in) :: a
    type(matrix3), intent(in) :: b
    type(matrix3) :: c
    !DIR$ FMA
    c%a11 = a%a11*b%a11 + a%a21*b%a21 + a%a31*b%a31
    c%a12 = a%a11*b%a12 + a%a21*b%a22 + a%a31*b%a32
    c%a13 = a%a11*b%a13 + a%a21*b%a23 + a%a31*b%a33
    c%a21 = a%a21*b%a11 + a%a22*b%a21 + a%a32*b%a31
    c%a22 = a%a21*b%a12 + a%a22*b%a22 + a%a32*b%a32
    c%a23 = a%a21*b%a13 + a%a22*b%a23 + a%a32*b%a33
    c%a31 = a%a31*b%a11 + a%a32*b%a21 + a%a33*b%a31
    c%a32 = a%a31*b%a12 + a%a32*b%a22 + a%a33*b%a32
    c%a33 = a%a31*b%a13 + a%a32*b%a23 + a%a33*b%a33    
    end function
    pure function smat3_product_smat3(a,b) result(c)
    type(symmatrix3), intent(in) :: a
    type(symmatrix3), intent(in) :: b
    type(matrix3) :: c
    !DIR$ FMA
    c%a11 = a%a11*b%a11 + a%a21*b%a21 + a%a31*b%a31
    c%a12 = a%a11*b%a21 + a%a21*b%a22 + a%a31*b%a32
    c%a13 = a%a11*b%a31 + a%a21*b%a32 + a%a31*b%a33
    c%a21 = a%a21*b%a11 + a%a22*b%a21 + a%a32*b%a31
    c%a22 = a%a21*b%a21 + a%a22*b%a22 + a%a32*b%a32
    c%a23 = a%a21*b%a31 + a%a22*b%a32 + a%a32*b%a33
    c%a31 = a%a31*b%a11 + a%a32*b%a21 + a%a33*b%a31
    c%a32 = a%a31*b%a21 + a%a32*b%a22 + a%a33*b%a32
    c%a33 = a%a31*b%a31 + a%a32*b%a32 + a%a33*b%a33    
    end function

    pure function vec3_cross_matrix(a) result(m)
    type(vector3), intent(in) :: a
    type(matrix3) :: m
    m%a11 = 0.0_wp
    m%a12 = -a%z
    m%a13 = a%y
    m%a21 = a%z
    m%a22 = 0.0_wp
    m%a23 = -a%x
    m%a31 = -a%y
    m%a32 = a%x
    m%a33 = 0.0_wp
    end function

    pure function vec3_outer(a,b) result(m)
    type(vector3), intent(in) :: a,b
    type(matrix3) :: m
    m%a11 = a%x*b%x
    m%a12 = a%x*b%y
    m%a13 = a%x*b%z
    m%a21 = a%y*b%x
    m%a22 = a%y*b%y
    m%a23 = a%y*b%z
    m%a31 = a%z*b%x
    m%a32 = a%z*b%y
    m%a33 = a%z*b%z
    end function

    pure function vec3_mmoi(v, negative) result(m)
    type(vector3), intent(in) :: v
    logical, intent(in), optional :: negative
    type(matrix3) :: m
    real(wp) :: xx,yy,zz,xy,zx,yz
    integer :: sign

    !tex: Inertia matrix
    ! $$-[\vec{v}\times][\vec{v}\times] = \begin{bmatrix}y^{2}+z^{2} & -xy & -zx\\
    !-xy & x^{2}+z^{2} & -yz\\
    !-zx & -yz & x^{2}+z^{2}
    !\end{bmatrix}$$

    if(present(negative) .and. negative) then
        sign = -1
    else
        sign  = 1
    end if
    !DIR$ FMA
    xx =  sign*v%x**2
    yy =  sign*v%y**2
    zz =  sign*v%z**2
    xy = -sign*v%x*v%y
    yz = -sign*v%y*v%z
    zx = -sign*v%z*v%x

    m%a11 = yy+zz
    m%a22 = xx+zz
    m%a33 = yy+xx
    m%a12 = xy
    m%a13 = zx
    m%a23 = yz
    m%a21 = xy
    m%a31 = zx
    m%a32 = yz

    end function

    pure function mat3_transpose(m) result(w)
    class(matrix3), intent(in) :: m
    type(matrix3) :: w
    w%a11 = m%a11
    w%a22 = m%a22
    w%a33 = m%a33
    w%a12 = m%a21
    w%a13 = m%a31
    w%a23 = m%a32
    w%a21 = m%a12
    w%a31 = m%a13
    w%a32 = m%a23
    end function

    pure function smat3_transpose(m) result(w)
    class(symmatrix3), intent(in) :: m
    type(symmatrix3) :: w
    w = m
    end function

    pure function mat3_get_column(m, col) result(v)
    class(matrix3), intent(in) :: m
    integer, intent(in) :: col
    type(vector3) :: v
    select case(col)
    case(1)
        v = vector3(m%a11, m%a21, m%a31)
    case(2)
        v = vector3(m%a12, m%a22, m%a32)
    case(3)
        v = vector3(m%a13, m%a23, m%a33)
    case default
        v = o_
    end select
    end function
    
    pure function mat3_get_row(m, row) result(v)
    class(matrix3), intent(in) :: m
    integer, intent(in) :: row
    type(vector3) :: v
    select case(row)
    case(1)
        v = vector3(m%a11, m%a12, m%a13)
    case(2)
        v = vector3(m%a21, m%a22, m%a23)
    case(3)
        v = vector3(m%a31, m%a32, m%a33)
    case default
        v = o_
    end select
    end function

    pure function smat3_get_column(m, col) result(v)
    class(symmatrix3), intent(in) :: m
    integer, intent(in) :: col
    type(vector3) :: v
    select case(col)
    case(1)
        v = vector3(m%a11, m%a21, m%a31)
    case(2)
        v = vector3(m%a21, m%a22, m%a32)
    case(3)
        v = vector3(m%a31, m%a32, m%a33)
    case default
        v = o_
    end select
    end function
    
    pure function smat3_get_row(m, row) result(v)
    class(symmatrix3), intent(in) :: m
    integer, intent(in) :: row
    type(vector3) :: v
    select case(row)
    case(1)
        v = vector3(m%a11, m%a21, m%a31)
    case(2)
        v = vector3(m%a21, m%a22, m%a32)
    case(3)
        v = vector3(m%a31, m%a32, m%a33)
    case default
        v = o_
    end select
    end function
    
    end module