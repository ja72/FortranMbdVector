!****************************************************************************
!
!  PROGRAM: FortranMbdVector
!
!  PURPOSE: Rigid Body Simulation using custom vector types
!
!****************************************************************************

    program FortranMbdVector
    use mod_vectors_dynamics
    implicit none

    ! Variables
    
    call mbd_solver(360*2520)

    end program FortranMbdVector

