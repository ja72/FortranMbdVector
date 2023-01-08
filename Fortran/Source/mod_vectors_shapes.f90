    module mod_vectors_shapes
    use mod_vectors_matrix3
    implicit none
    
    ! Supported geometry tags
    enum, bind(c)
        enumerator :: sphere
        enumerator :: cylinder
    end enum

    ! Generic geometry for body
    type, abstract :: shapes
        ! Contain mass and mass moment of inertia principal values
        real(wp) :: volume
        type(vector3) :: center
        real(wp) :: vmmoi(3)
    contains
        procedure(nearest_point_function), deferred, pass :: nearest_point        
    end type

    ! Geometry of a sphere
    type, extends(shapes) :: sphere_shape
        real(wp) :: radius
    contains
        procedure :: nearest_point => nearest_point_sphere
    end type
    
    ! Geometry of cylinder
    type, extends(shapes) :: cylinder_shape
        real(wp) :: length, radius
    contains
        procedure :: nearest_point => nearest_point_cylinder
    end type

    interface shapes
        procedure :: new_shape
    end interface

    abstract interface
        pure function nearest_point_function(shape, direction) result(point)
        ! Finds the point of the surface of the body closest to the plane with specified
        ! normal direction. The point is given relative to the center of the body in
        ! world coordinates.
        import
        class(shapes), intent(in) :: shape
        type(vector3), intent(in) :: direction
        type(vector3) :: point
        end function
    end interface
    
    
    contains
    
    
    pure function new_shape(shape_type, dims) result(s)
    ! Constructor for body geometry based on tag and dimensions
    integer, intent(in) :: shape_type
    real(wp), intent(in) :: dims(:)
    class(shapes), allocatable :: s
    real(wp) :: V, r, l, vmmoi(3)
        select case (shape_type)
        case (sphere)
            r = dims(1)
            !tex: Volume of sphere $V = \frac{4}{3} \pi r^3$
            V = 4*pi*r**3/3
            !tex: MMMOI of sphere $I_{\rm all}=\frac{2}{5} m r^2$
            vmmoi = [ 2*r**2/5, 2*r**2/5, 2*r**2/5 ]
            s = sphere_shape(V, vector3(0.0_wp,0.0_wp,0.0_wp), vmmoi, r)
        case (cylinder)                
            l = dims(1)
            r = dims(2)
            !tex: Volume of cylinder $V = \ell \pi r^2$
            V = l*pi*r**2
            !tex: MMMOI of cylinder $I_{\rm zz}=\frac{m}{2} r^2$, $I_{\rm xx}=\frac{m}{4} r^2 + \frac{m}{12} \ell^2$
            vmmoi = [l**2/12 + r**2/4, l**2/12 + r**2/4, r**2/2]
            s = cylinder_shape(V, vector3(0.0_wp,0.0_wp,0.0_wp), vmmoi, l, r)
        end select
    end function

    pure function nearest_point_sphere(shape, direction) result(point)
    ! Find the point of the surface of the body closest to the plane with specified
    ! normal direction. The point is given relative to the center of the body in
    ! world coordinates.
    ! If for example the horizontal plane has normal [0,0,1] then
    ! the nearest point on the sphere to the plane is [0,0,-r] 
    class(sphere_shape), intent(in) :: shape
    type(vector3), intent(in) :: direction
    type(vector3) :: point
        !tex: Point $\vec{\rm pos} = -r\, \vec{e}$
        point = -direction * shape%radius
    end function
    
    pure function nearest_point_cylinder(shape, direction) result(point)
    ! Find the point of the surface of the body closest to the plane with specified
    ! normal direction. The point is given relative to the center of the body in
    ! world coordinates.
    ! For a cylinder the nearest point has to be on the circle the defines the
    ! end caps. The orientation angle of the local x-axis towards the plane is 
    ! calculated first. That locates the point around the cylinder and the end
    ! cap is chosen by the dot product of the local z-axis to the plane.
    class(cylinder_shape), intent(in) :: shape
    type(vector3), intent(in) :: direction
    type(vector3) :: point
    real(wp) :: d,l,phi,alp,x,y
    
        !tex: Point $\vec{\rm pos} = \pmatrix{
        !  \frac{d}{2} \cos \varphi \\ 
        !  \frac{d}{2} \sin \varphi \\ 
        !  \pm \frac{\ell}{2} }$
        d = shape%radius*2
        l = shape%length
        alp = -sign(1.0, direction%z)
        x = direction%y
        y = direction%z
        if( abs(x)<=tiny .and. abs(y)<=tiny) then
            point = vector3(0.0_wp, 0.0_wp, (l/2)*alp)
        else
            phi = atan2(y, x) - pi
            point = vector3((d/2)*cos(phi), (d/2)*sin(phi), (l/2)*alp)
        end if
        
    end function
    
    
    end module