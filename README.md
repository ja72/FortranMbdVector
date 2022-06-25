# FortranMbdVector

Fortran application for the dynamic simulation a single rigid body. Some of the highlights of the code are:

 - Momentum based state vector.
 - Handles a single contact with the floor, including coefficient of restitution and friction.
 - Custom vector/matrix/quaternion types
 - Speed is about 2M simulation steps per second.
 - Mass properties derived from shape definition (also used in contacts)

Application is also ported to CSharp for comparison of speed and accuracy.

### Screenshot

![image](https://user-images.githubusercontent.com/22509289/175482233-1c1a6f26-f70d-427c-875c-7e36c8330258.png)

### Fortraan User Types

```fortran

    ! Standard vector
    type :: vector3
        real(wp) :: x,y,z
    contains
        procedure :: to_array => vec3_to_array
    end type
    
    ! Quaternion for rotation encoding
    type :: quaternion
        type(vector3) :: vector
        real(wp) :: scalar
    contains
        procedure :: to_array => q_to_array
    end type
    
    ! Standard 3×3 matrix
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
    
    ! Generic geometry for body
    type, abstract :: shapes
        ! Contain mass and mass moment of inertia principal values
        real(wp) :: volume
        type(vector3) :: center
        real(wp) :: vmmoi(3)
    contains
        procedure(nearest_point_function), deferred, pass :: nearest_point        
    end type
    
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
    
```

