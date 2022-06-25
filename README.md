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

### Fortran User Types

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

## Rigid Body Theory

### Kinematics

I am using the momentum vectors to keep track of the motion of the body because they are easier to integrate over time. As summed at the center of mass, based on the orientation quaterion $q$ the following equations estimate the mass moment of inertia tensor

$$
\begin{aligned}{\rm R} & ={\rm rot}\left(\mathscr{q}\right)\\
{\bf I}_{c} & ={\rm R}{\bf I}_{body}{\rm R}^{\intercal}
\end{aligned}
$$

and the momentum vectors of a rigid body are

$$
\begin{aligned}\boldsymbol{p} & =m\boldsymbol{v}_{c}\\
\boldsymbol{L}_{c} & ={\bf I}_{c}\boldsymbol{\omega}
\end{aligned}
$$

To extract the motion vectors from the above momentum use reverse of the above

$$
\begin{aligned}\boldsymbol{v}_{c} & =\frac{1}{m}\boldsymbol{p}\\
\boldsymbol{\omega} & ={\bf I}_{c}^{-1}\boldsymbol{L}_{c}
\end{aligned}
$$

As a generalization the refrence point **b** being tracked is not the center of mass. At each time step a vector $\boldsymbol{c}$ is calculating the holds the location of the center of mass relative to the reference point. Now the following transformations are needed 

$$
\begin{aligned}\boldsymbol{v}_{c} & =\boldsymbol{v}_{b}+\boldsymbol{\omega}\times\boldsymbol{c}\\
\boldsymbol{L}_{b} & =\boldsymbol{L}_{c}+\boldsymbol{c}\times\boldsymbol{p}
\end{aligned}
$$

And to extract the motion vectors use

$$
\begin{aligned}\boldsymbol{v}_{b} & =\frac{1}{m}\boldsymbol{p}+\boldsymbol{c}\times\boldsymbol{\omega}\\
\boldsymbol{\omega} & ={\bf I}_{c}^{-1}\left(\boldsymbol{L}_{b}-\boldsymbol{c}\times\boldsymbol{p}\right)
\end{aligned}
$$

### Dynamics

Newton's 2nd law of motion provides the link between forces and momenum as expressed on the center of mass

$$\begin{aligned}\boldsymbol{F} & =\frac{{\rm d}}{{\rm d}t}\boldsymbol{p}\\
\boldsymbol{\tau}_{c} & =\frac{{\rm d}}{{\rm d}t}\boldsymbol{L}_{c}
\end{aligned}$$

but when expressed at a different reference point the above vecomes a bit more complex

$$
\begin{aligned}\boldsymbol{F} & =\frac{{\rm d}}{{\rm d}t}\boldsymbol{p}\\
\boldsymbol{\tau}_{b} & =\frac{{\rm d}}{{\rm d}t}\left(\boldsymbol{L}_{b}\right)+\boldsymbol{v}_{b}\times\boldsymbol{p}
\end{aligned}
$$

### Body State Vector

The combined position, orientation and momentum vectors define the body state vector $Y$

$$
Y=\begin{Bmatrix}\boldsymbol{r}_{b}\\
\mathscr{q}\\
\boldsymbol{p}\\
\boldsymbol{L}_{b}
\end{Bmatrix}$$

And from the above dynamics at each time step the following is integrated to produce the next body state $Y \rightarrow Y + h \frac{{\rm d}}{{\rm d}t}Y$

$$\frac{{\rm d}}{{\rm d}t}Y=\begin{Bmatrix}\dot{\boldsymbol{r}}_{b}\\
\dot{\mathscr{q}}\\
\dot{\boldsymbol{p}}\\
\dot{\boldsymbol{L}}_{b}
\end{Bmatrix}=\begin{Bmatrix}\frac{1}{m}\boldsymbol{p}+\boldsymbol{c}\times\boldsymbol{\omega}\\
\tfrac{1}{2}\boldsymbol{\omega}\otimes\mathscr{q}\\
\boldsymbol{F}\\
\boldsymbol{\tau}_{b}-\boldsymbol{v}_{b}\times\boldsymbol{p}
\end{Bmatrix}$$

The actual integrator in the code is a simple Runge-Kutta 4 method implemented as a type bound procedure to the `simulation` object

```fortran
    pure function sim_integrate(sim, h) result(Y_next)
    ! Implement RK4 integrator for rigid bodies
    class(simulation), intent(in) :: sim
    real(wp), intent(in) :: h
    real(wp), dimension(state_size) :: Y_next, Y1, Y2, Y3
    real(wp) :: k0(state_size),k1(state_size),k2(state_size),k3(state_size)
    
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
```
