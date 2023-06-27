!==========================================================================================!
!                           thermal convection using GDS-WCSPH                             !
!------------------------------------------------------------------------------------------!
!                          Copyright by Kensuke Shobuzako (2022)                           !
!==========================================================================================!

module INPUT

    implicit none

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    !                                     INPUT below                                      !
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    !================================== system setting ====================================!
    ! NAME
    character(99), parameter :: file_name = 'R100'       ! file I.D.
    ! RSST_VIM system
    integer, parameter :: RSST_VIM        = 0            ! 0(air), 1(mantle)
    ! the number of threads(cores)
    integer, parameter :: num_threads     = 8            ! the number of threads
    ! STEPS
    integer, parameter :: total_step      = 500          ! total steps
    integer, parameter :: write_step      = 10           ! writing step
    real(8), parameter :: coe_step        = 0.3d0        ! coefficient of time step 
    ! SYSTEM INFO
    real(8), parameter :: Ra_0            = 1.0d+4       ! Ra_0 number
    integer, parameter :: num             = 50           ! the number of particles in an edge
    ! RELAXATION PARAMETER
    real(8), parameter :: zeta            = 10800.d0     ! relaxation parameter (RSST)
    real(8), parameter :: xi              = 1.0d0        ! relaxation parameter (VIM)
    ! the number of analysis point (default:50)
    integer, parameter :: num_MP_side     = 100          ! the number of MP in an edge
    !================================= kernel function ====================================!
    integer, parameter :: kernel_switch   = 2            ! 0(Quintic), 1(WenC4), 2(WenC6)
    real(8), parameter :: coe_smooth      = 2.3d0        ! coefficient of smoothing length
    !================================ boudary condition ===================================!
    integer, parameter :: wall_v          = 0            ! 0(no-slip), 1(free-slip) 
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
    !                                     INPUT end                                        !
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!


    !============================== calculation parameter =================================!
    real(8), parameter :: g        = -10.0d0                ! gravity [m s-2]
    ! for air convection (Pr < Ra)
    real(8), parameter :: length_R = 1.0d-1                ! an edge length [m]
    real(8), parameter :: T_top_R  = 298.0d0               ! temperature of the top [K]
    real(8), parameter :: rho_0_R  = 1.0d0                 ! reference density [kg m-3]
    real(8), parameter :: k_h_R    = 1.0d-2                ! heat conductivity [W m-1 K-1]
    real(8), parameter :: eta_R    = 7.1d-6                ! viscosity [Pa s]
    real(8), parameter :: alpha_R  = 7.1d-3                ! thermal expansion [K-1]
    real(8), parameter :: c_p_R    = 1.0d+3                ! specific heat [J kg-1 K-1]
    real(8), parameter :: K_0_R    = 1.4d+5                ! bulk modulus [Pa]
    !======================================================================================!
    ! for mantle convection (Pr > Ra); Blankenbach (1989)
    real(8), parameter :: length_V = 1.0d+6                ! an edge length [m]
    real(8), parameter :: T_top_V  = 0.0d0                 ! temperature of the top [K]
    real(8), parameter :: T_bot_V  = 1.0d+3                ! temperature of the bottom [K]
    real(8), parameter :: rho_0_V  = 4.0d+3                ! reference density [kg m-3]
    real(8), parameter :: k_h_V    = 5.0d0                 ! heat conductivity [W m-1 K-1]
    real(8), parameter :: alpha_V  = 2.5d-5                ! thermal expansion [K-1]
    real(8), parameter :: c_p_V    = 1250.0d0              ! specific heat [J kg-1 K-1]
    real(8), parameter :: K_0_V    = 100.0d+9              ! bulk modulus [Pa]
    !======================================================================================!

end module INPUT