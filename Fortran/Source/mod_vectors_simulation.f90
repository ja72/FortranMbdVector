    module mod_vectors_simulation
    use mod_vectors_contact
    implicit none
    
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