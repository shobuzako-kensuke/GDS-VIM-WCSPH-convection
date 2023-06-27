!==========================================================================================!
!                           thermal convection using GDS-WCSPH                             !
!------------------------------------------------------------------------------------------!
!                          Copyright by Kensuke Shobuzako (2022)                           !
!==========================================================================================!

module set_time_step

    use INPUT
    use GLB_SETTING

    real(8) :: tmp_CFL, tmp_mo

contains
    subroutine set_delta_t
        ! fixed time step
        if (t == 1) then
            time_step = 0.01d0 * min(r_t_CFL, r_t_VNM, t_VNT)
        else
            time_step = min(r_t_CFL, r_t_VNM, t_VNT)
        endif
        pass_time = pass_time + time_step
    end subroutine set_delta_t

end module set_time_step