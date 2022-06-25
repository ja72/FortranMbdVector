    module mod_constants
    use, intrinsic :: iso_fortran_env, only: int64, real32, real64, real128
    use ifcore
    implicit none
    
    ! Working pression
    integer, parameter :: qw = int64, wp = real64

    ! Common constants
    real(wp), parameter :: pi = 3.1415926535897932384626433832795_wp
    real(wp), parameter :: deg = pi/180, mm=1/1000.0_wp, inch=0.0254_wp;
    real(wp), parameter :: ulp = 1.0_wp/2251799813685248
    real(wp), parameter :: tiny = 64*ulp, small = 1.0_wp/8388608

    end module