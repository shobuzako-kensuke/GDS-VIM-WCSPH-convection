!==========================================================================================!
!                           thermal convection using GDS-WCSPH                             !
!------------------------------------------------------------------------------------------!
!                          Copyright by Kensuke Shobuzako (2022)                           !
!==========================================================================================!

module analysis
    use INPUT
    use GLB_SETTING
    use GLB_ARG
    use INITIAL_SETTING
    use NNPS

    implicit none
    real(8), save :: V_rms_ave

contains
    subroutine analysis_exe
        tmp = 0.0d0
        do i = 1, total_num_system
            if (SP_kind(i) == 0) then
                ! cal V_rms
                tmp = tmp + (SP_uv(i, 1)**2 + SP_uv(i, 2)**2)
            endif
        enddo
        V_rms_ave = dsqrt(tmp / dble(total_num))
    end subroutine analysis_exe
end module analysis