!==========================================================================================!
!                           thermal convection using GDS-WCSPH                             !
!------------------------------------------------------------------------------------------!
!                          Copyright by Kensuke Shobuzako (2022)                           !
!==========================================================================================!

module cal_KERNEL

    use INPUT
    use GLB_SETTING

    implicit none
    real(8), save :: r_ij, q_ij, W_ij, W_ij_ave, dW_ij

contains
    !---------------------
    ! W_ave for PS scheme
    !---------------------
    subroutine cal_W_ave
        r_ij = delta_x
        q_ij = r_ij / smooth_len
        ! Quintic Spline
        if (kernel_switch == 0) then
            if ((0.0d0 <= q_ij) .and. (q_ij < 1.0d0)) then
                W_ij  = alpha_d * (            (3.0d0 - q_ij)**5 &
                                    - 6.0d0  * (2.0d0 - q_ij)**5 &
                                    + 15.0d0 * (1.0d0 - q_ij)**5 )

            elseif ((1.0d0 <= q_ij) .and. (q_ij < 2.0d0)) then
                W_ij  = alpha_d * (           (3.0d0 - q_ij)**5 &
                                    - 6.0d0 * (2.0d0 - q_ij)**5 )

            elseif ( ( 2.0d0 <= q_ij ) .and. ( q_ij <= 3.0d0 ) ) then
                W_ij  = alpha_d * (3.0d0 - q_ij)**5

            else
                W_ij  = 0.0d0
            end if
        
        ! Wendland C4
        elseif (kernel_switch == 1) then
            if (q_ij <= 2.0d0) then
                W_ij = alpha_d * ((1.0d0 - 0.5d0*q_ij)**6 &
                        * (1.0d0 + 3.0d0*q_ij + 35.0d0/12.0d0*q_ij**2))
            else
                W_ij = 0.0d0
            end if    

        ! Wendland C6
        elseif (kernel_switch == 2) then
            if (q_ij <= 2.0d0) then
                W_ij = alpha_d * ((1.0d0 - 0.5d0*q_ij)**8 &
                        * (1.0d0 + 4.0d0*q_ij + 25.0d0/4.0d0*q_ij**2 + 4.0d0*q_ij**3))
            else
                W_ij = 0.0d0
            end if
        end if

        W_ij_ave = W_ij
    end subroutine cal_W_ave

    !---------------------
    ! cal W_ij and dW_ij
    !---------------------
    subroutine cal_W_ij_and_dW_ij(q_ij, W_ij, dW_ij)
        implicit none
        real(8), intent(in)  :: q_ij
        real(8), intent(out) :: W_ij, dW_ij

        ! Quintic Spline
        if (kernel_switch == 0) then
            if ((0.0d0 <= q_ij) .and. (q_ij < 1.0d0)) then
                W_ij  = alpha_d * (            (3.0d0 - q_ij)**5 &
                                    - 6.0d0  * (2.0d0 - q_ij)**5 &
                                    + 15.0d0 * (1.0d0 - q_ij)**5 )
                dW_ij = - (5.0d0 / smooth_len) * alpha_d &
                                * (            (3.0d0 - q_ij)**4 &
                                    - 6.0d0  * (2.0d0 - q_ij)**4 &
                                    + 15.0d0 * (1.0d0 - q_ij)**4 )

            elseif ((1.0d0 <= q_ij) .and. (q_ij < 2.0d0)) then
                W_ij  = alpha_d * (           (3.0d0 - q_ij)**5 &
                                    - 6.0d0 * (2.0d0 - q_ij)**5 )
                dW_ij = - (5.0d0 / smooth_len) * alpha_d &
                               * (            (3.0d0 - q_ij)**4 &
                                   - 6.0d0  * (2.0d0 - q_ij)**4 )

            elseif ((2.0d0 <= q_ij) .and. (q_ij <= 3.0d0)) then
                W_ij  = alpha_d * (3.0d0 - q_ij)**5
                dW_ij = - (5.0d0 / smooth_len) * alpha_d * (3.0d0 - q_ij)**4

            else
                W_ij  = 0.0d0
                dW_ij = 0.0d0
            end if

        ! Wendland C4
        elseif (kernel_switch == 1) then
            if (q_ij <= 2.0d0) then
                W_ij = alpha_d * ((1.0d0 - 0.5d0*q_ij)**6 &
                        * (1.0d0 + 3.0d0*q_ij + 35.0d0/12.0d0*q_ij**2))
                dW_ij = 3.0d0 * alpha_d / smooth_len &
                        * ( - (1.0d0 - 0.5d0*q_ij)**5 * (1.0d0 + 3.0d0*q_ij + 35.0d0/12.0d0*q_ij**2) &
                            + (1.0d0 - 0.5d0*q_ij)**6 * (1.0d0 + 35.0d0/18.0d0*q_ij)                 )
            else
                W_ij  = 0.0d0
                dW_ij = 0.0d0
            end if

        ! Wendland C6
        elseif (kernel_switch == 2) then
            if (q_ij <= 2.0d0) then
                W_ij = alpha_d * ((1.0d0 - 0.5d0*q_ij)**8 &
                        * (1.0d0 + 4.0d0*q_ij + 25.0d0/4.0d0*q_ij**2 + 4.0d0*q_ij**3))
                dW_ij = 4.0d0 * alpha_d / smooth_len &
                        * ( - (1.0d0 - 0.5d0*q_ij)**7 * (1.0d0 + 4.0d0*q_ij + 25.0d0/4.0d0*q_ij**2 + 4.0d0*q_ij**3) &
                            + (1.0d0 - 0.5d0*q_ij)**8 * (1.0d0 + 25.0d0/8.0d0*q_ij + 3.0d0*q_ij**2)                 )
            else
                W_ij  = 0.0d0
                dW_ij = 0.0d0
            end if
        end if
    end subroutine cal_W_ij_and_dW_ij

    !---------------------
    ! cal W_ij
    !---------------------
    subroutine cal_W_ij(q_ij, W_ij)
        implicit none
        real(8), intent(in)  :: q_ij
        real(8), intent(out) :: W_ij

        ! Quintic Spline
        if (kernel_switch == 0) then
            if ((0.0d0 <= q_ij) .and. (q_ij < 1.0d0)) then
                W_ij  = alpha_d * (            (3.0d0 - q_ij)**5 &
                                    - 6.0d0  * (2.0d0 - q_ij)**5 &
                                    + 15.0d0 * (1.0d0 - q_ij)**5 )
                
            elseif ((1.0d0 <= q_ij) .and. (q_ij < 2.0d0)) then
                W_ij  = alpha_d * (           (3.0d0 - q_ij)**5 &
                                    - 6.0d0 * (2.0d0 - q_ij)**5 )
                
            elseif ((2.0d0 <= q_ij) .and. (q_ij <= 3.0d0)) then
                W_ij  = alpha_d * (3.0d0 - q_ij)**5
                
            else
                W_ij  = 0.0d0
            end if

        ! Wendland C4
        elseif (kernel_switch == 1) then
            if (q_ij <= 2.0d0) then
                W_ij = alpha_d * ((1.0d0 - 0.5d0*q_ij)**6 &
                        * (1.0d0 + 3.0d0*q_ij + 35.0d0/12.0d0*q_ij**2))
            else
                W_ij  = 0.0d0
            end if

        ! Wendland C6
        elseif (kernel_switch == 2) then
            if (q_ij <= 2.0d0) then
                W_ij = alpha_d * ((1.0d0 - 0.5d0*q_ij)**8 &
                        * (1.0d0 + 4.0d0*q_ij + 25.0d0/4.0d0*q_ij**2 + 4.0d0*q_ij**3))
            else
                W_ij  = 0.0d0
            end if
        end if
    end subroutine cal_W_ij

    !---------------------
    ! cal dW_ij
    !---------------------
    subroutine cal_dW_ij(q_ij, dW_ij)
        implicit none
        real(8), intent(in)  :: q_ij
        real(8), intent(out) :: dW_ij

        ! Quintic Spline
        if (kernel_switch == 0) then
            if ((0.0d0 <= q_ij) .and. (q_ij < 1.0d0)) then
                dW_ij = - (5.0d0 / smooth_len) * alpha_d &
                                * (            (3.0d0 - q_ij)**4 &
                                    - 6.0d0  * (2.0d0 - q_ij)**4 &
                                    + 15.0d0 * (1.0d0 - q_ij)**4 )

            elseif ((1.0d0 <= q_ij) .and. (q_ij < 2.0d0)) then
                dW_ij = - (5.0d0 / smooth_len) * alpha_d &
                               * (            (3.0d0 - q_ij)**4 &
                                   - 6.0d0  * (2.0d0 - q_ij)**4 )

            elseif ((2.0d0 <= q_ij) .and. (q_ij <= 3.0d0)) then
                dW_ij = - (5.0d0 / smooth_len) * alpha_d * (3.0d0 - q_ij)**4

            else
                dW_ij = 0.0d0
            end if

        ! Wendland C4
        elseif (kernel_switch == 1) then
            if (q_ij <= 2.0d0) then
                dW_ij = 3.0d0 * alpha_d / smooth_len &
                        * ( - (1.0d0 - 0.5d0*q_ij)**5 * (1.0d0 + 3.0d0*q_ij + 35.0d0/12.0d0*q_ij**2) &
                            + (1.0d0 - 0.5d0*q_ij)**6 * (1.0d0 + 35.0d0/18.0d0*q_ij)                 )
            else
                dW_ij = 0.0d0
            end if

        ! Wendland C6
        elseif (kernel_switch == 2) then
            if (q_ij <= 2.0d0) then
                dW_ij = 4.0d0 * alpha_d / smooth_len &
                        * ( - (1.0d0 - 0.5d0*q_ij)**7 * (1.0d0 + 4.0d0*q_ij + 25.0d0/4.0d0*q_ij**2 + 4.0d0*q_ij**3) &
                            + (1.0d0 - 0.5d0*q_ij)**8 * (1.0d0 + 25.0d0/8.0d0*q_ij + 3.0d0*q_ij**2)                 )
            else
                dW_ij = 0.0d0
            end if
        end if
    end subroutine cal_dW_ij
    
end module cal_KERNEL