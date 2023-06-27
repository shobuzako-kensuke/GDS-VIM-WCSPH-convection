!==========================================================================================!
!                           thermal convection using GDS-WCSPH                             !
!------------------------------------------------------------------------------------------!
!                          Copyright by Kensuke Shobuzako (2022)                           !
!==========================================================================================!

module cal_EOS

    use INPUT
    use GLB_SETTING

    implicit none

    real(8) :: rho_S, pre
    real(8) :: tem, rho_G

contains

    !---------------------
    ! artificial EOS
    !---------------------
    subroutine cal_artificial_EOS(rho_S, pre)
        real(8), intent(in)  :: rho_S
        real(8), intent(out) :: pre

        pre = (K_0 / rho_0) / (zeta*zeta) * (rho_S - rho_0)

    end subroutine cal_artificial_EOS
    
    !---------------------
    ! real EOS
    !---------------------
    subroutine cal_real_EOS(pre, tem, rho_G)
        real(8), intent(in)  :: pre, tem
        real(8), intent(out) :: rho_G

        rho_G = rho_0 * (- alpha * (tem - T_0) + pre / K_0)
    end subroutine cal_real_EOS

end module cal_EOS