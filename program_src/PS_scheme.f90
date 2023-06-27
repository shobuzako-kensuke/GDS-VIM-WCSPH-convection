!==========================================================================================!
!                           thermal convection using GDS-WCSPH                             !
!------------------------------------------------------------------------------------------!
!                          Copyright by Kensuke Shobuzako (2022)                           !
!==========================================================================================!

module PS_scheme
    
    use INPUT
    use GLB_SETTING
    use GLB_ARG
    use cal_EOS
    use INITIAL_SETTING
    use NNPS
    use cal_KERNEL

    implicit none
    integer :: my_cell, your_cell
    real(8) :: u_max_tmp, v_max_tmp, V_max, x_ij, y_ij, dW_ij_x, dW_ij_y, u_ij, v_ij
    real(8) :: coe_PS, KGC_dW_ij_x, KGC_dW_ij_y

contains
    subroutine PS_scheme_exe
        PS_shift_xy(:,:)  = 0.0d0
        PS_TE_uv(:,:)     = 0.0d0
        PS_TE_rho_S(:,:)  = 0.0d0
        PS_TE_pre(:,:)    = 0.0d0
        PS_TE_tem(:,:)    = 0.0d0
        PS_TE_sm_rho(:,:) = 0.0d0
        PS_KGC(:,:)       = 0.0d0

        ! set PS coefficient
        u_max_tmp = maxval(SP_uv(:, 1))
        v_max_tmp = maxval(SP_uv(:, 2))
        V_max = max(u_max_tmp, v_max_tmp)

        coe_PS = 0.5d0 * V_max * time_step * smooth_len

        !$ call omp_set_dynamic(.false.)
        !$ call omp_set_num_threads(num_threads)
        !$OMP parallel default(none) &
        !$OMP private(i, j, k, me, you, my_cell, your_cell) &
        !$OMP private(x_ij, y_ij, r_ij, q_ij, W_ij, dW_ij) &
        !$OMP private(tmp, dW_ij_x, dW_ij_y, u_ij, v_ij) &
        !$OMP private(KGC_dW_ij_x, KGC_dW_ij_y) &
        !$OMP shared(inner_cell_start, inner_cell_end, total_num_system) &
        !$OMP shared(NNPS_arg, SP_kind, NNPS_9, NNPS_max) &
        !$OMP shared(SP_xy, SP_uv, SP_rho_S, SP_pre, SP_tem, SP_sm_rho, SP_rho_G) &
        !$OMP shared(smooth_len, coe_PS, SP_mass) &
        !$OMP shared(PS_shift_xy, PS_TE_uv, PS_TE_rho_S, PS_TE_pre, PS_TE_tem, PS_TE_sm_rho) &
        !$OMP shared(PS_KGC, W_ij_ave, NNPS_row, time_step, V_max) &
        !$OMP shared(wall_width, length) 

        !-------------------------
        ! cal KGC
        !-------------------------
        !$OMP do
        ! my cell number loop
        do my_cell = inner_cell_start, inner_cell_end

            ! my number loop
            do i = 1, NNPS_row(my_cell)-1
                me = NNPS_arg(i, my_cell) ! my number

                ! if I am inner_particle, exe
                if (SP_kind(me) == 0) then

                    
                    ! your cell loop
                    do j = 1, 9
                        your_cell = NNPS_9(j, my_cell)

                        ! your number loop
                        do k = 1, NNPS_row(your_cell)-1
                            you = NNPS_arg(k, your_cell)

                            ! if you are not me, exe
                            if (me /= you) then
                                !-----------------------
                                ! cal kernel
                                !-----------------------
                                x_ij = SP_xy(me, 1) - SP_xy(you, 1)
                                y_ij = SP_xy(me, 2) - SP_xy(you, 2)
                                r_ij = dsqrt(x_ij*x_ij + y_ij*y_ij)
                                q_ij = r_ij / smooth_len
                                call cal_dW_ij(q_ij, dW_ij)

                                !-----------------------
                                ! cal KGC
                                !-----------------------
                                tmp = SP_mass / SP_rho_S(you)
                                dW_ij_x = x_ij / r_ij * dW_ij
                                dW_ij_y = y_ij / r_ij * dW_ij
                                
                                PS_KGC(me, 1) = PS_KGC(me, 1) + tmp * (-x_ij) * dW_ij_x
                                PS_KGC(me, 2) = PS_KGC(me, 2) + tmp * (-y_ij) * dW_ij_x
                                PS_KGC(me, 3) = PS_KGC(me, 3) + tmp * (-x_ij) * dW_ij_y
                                PS_KGC(me, 4) = PS_KGC(me, 4) + tmp * (-y_ij) * dW_ij_y
                            endif
                        enddo
                    enddo
                endif
            enddo
        enddo
        !$OMP enddo
        !$OMP barrier

        !-----------------------
        ! cal inverse KGC arg
        !-----------------------
        !$OMP do
        do i = 1, total_num_system
            if (SP_kind(i) == 0) then
                tmp = 1.0d0 / (PS_KGC(i, 1)*PS_KGC(i, 4) - PS_KGC(i, 2)*PS_KGC(i, 3))
                x_ij = PS_KGC(i, 1)
                PS_KGC(i, 1) =   tmp * PS_KGC(i, 4)
                PS_KGC(i, 2) = - tmp * PS_KGC(i, 2)
                PS_KGC(i, 3) = - tmp * PS_KGC(i, 3)
                PS_KGC(i, 4) =   tmp * x_ij
            endif
        enddo
        !$OMP enddo
        !$OMP barrier

        !-------------------------
        ! cal PS & (nabla phi)_i
        !-------------------------
        !$OMP do
        ! my cell number loop
        do my_cell = inner_cell_start, inner_cell_end

            ! my number loop
            do i = 1, NNPS_row(my_cell)-1
                me = NNPS_arg(i, my_cell) ! my number

                ! if I am inner_particle, exe
                if (SP_kind(me) == 0) then

                    ! your cell loop
                    do j = 1, 9
                        your_cell = NNPS_9(j, my_cell)

                        ! your number loop
                        do k = 1, NNPS_row(your_cell)-1
                            you = NNPS_arg(k, your_cell)

                            ! if you are not me, exe
                            if (me /= you) then
                                !-----------------------
                                ! cal kernel
                                !-----------------------
                                x_ij = SP_xy(me, 1) - SP_xy(you, 1)
                                y_ij = SP_xy(me, 2) - SP_xy(you, 2)
                                r_ij = dsqrt(x_ij*x_ij + y_ij*y_ij)
                                q_ij = r_ij / smooth_len
                                call cal_W_ij_and_dW_ij(q_ij, W_ij, dW_ij)

                                !-----------------------
                                ! cal right-hands
                                !-----------------------
                                tmp = SP_mass / SP_rho_S(you)
                                dW_ij_x = x_ij / r_ij * dW_ij
                                dW_ij_y = y_ij / r_ij * dW_ij
                                KGC_dW_ij_x = PS_KGC(me, 1)*dW_ij_x + PS_KGC(me, 2)*dW_ij_y
                                KGC_dW_ij_y = PS_KGC(me, 3)*dW_ij_x + PS_KGC(me, 4)*dW_ij_y

                                u_ij = SP_uv(me, 1) - SP_uv(you, 1)
                                v_ij = SP_uv(me, 2) - SP_uv(you, 2)

                                ! cal_shift
                                PS_shift_xy(me, 1) = PS_shift_xy(me, 1) - coe_PS * tmp * KGC_dW_ij_x
                                PS_shift_xy(me, 2) = PS_shift_xy(me, 2) - coe_PS * tmp * KGC_dW_ij_y

                                ! Taylor expansion
                                PS_TE_uv(me, 1)    = PS_TE_uv(me, 1)    + tmp * (SP_uv(you, 1) - SP_uv(me, 1)) * KGC_dW_ij_x
                                PS_TE_uv(me, 2)    = PS_TE_uv(me, 2)    + tmp * (SP_uv(you, 1) - SP_uv(me, 1)) * KGC_dW_ij_y
                                PS_TE_uv(me, 3)    = PS_TE_uv(me, 3)    + tmp * (SP_uv(you, 2) - SP_uv(me, 2)) * KGC_dW_ij_x
                                PS_TE_uv(me, 4)    = PS_TE_uv(me, 4)    + tmp * (SP_uv(you, 2) - SP_uv(me, 2)) * KGC_dW_ij_y
                                PS_TE_rho_S(me, 1) = PS_TE_rho_S(me, 1) + tmp * (SP_rho_S(you) - SP_rho_S(me)) * KGC_dW_ij_x
                                PS_TE_rho_S(me, 2) = PS_TE_rho_S(me, 2) + tmp * (SP_rho_S(you) - SP_rho_S(me)) * KGC_dW_ij_y
                                PS_TE_tem(me, 1)   = PS_TE_tem(me, 1)   + tmp * (SP_tem(you)   - SP_tem(me))   * KGC_dW_ij_x
                                PS_TE_tem(me, 2)   = PS_TE_tem(me, 2)   + tmp * (SP_tem(you)   - SP_tem(me))   * KGC_dW_ij_y
                                ! PS_TE_pre(me, 1) = PS_TE_pre(me, 1) + tmp * (SP_pre(you) - SP_pre(me)) * KGC_dW_ij_x
                                ! PS_TE_pre(me, 2) = PS_TE_pre(me, 2) + tmp * (SP_pre(you) - SP_pre(me)) * KGC_dW_ij_y
                                ! PS_TE_sm_rho(me, 1) = PS_TE_sm_rho(me, 1) + tmp * (SP_sm_rho(you) - SP_sm_rho(me)) * KGC_dW_ij_x
                                ! PS_TE_sm_rho(me, 2) = PS_TE_sm_rho(me, 2) + tmp * (SP_sm_rho(you) - SP_sm_rho(me)) * KGC_dW_ij_y
                            endif
                        enddo
                    enddo
                endif
            enddo
        enddo
        !$OMP enddo
        !$OMP barrier

        !-----------------------
        ! shifting & renewal
        !-----------------------
        !$OMP do
        do i = 1, total_num_system
            if (SP_kind(i) == 0) then
                ! shift
                tmp = dsqrt(PS_shift_xy(i, 1)**2 + PS_shift_xy(i, 2)**2)

                if (abs(PS_shift_xy(i, 1)) >= 0.1d0*V_max*time_step) then
                    PS_shift_xy(i, 1) = PS_shift_xy(i, 1) / tmp * (0.1d0*V_max*time_step)
                    SP_xy(i, 1) = SP_xy(i, 1) + PS_shift_xy(i, 1)
                else
                    SP_xy(i, 1) = SP_xy(i, 1) + PS_shift_xy(i, 1)
                endif
    
                if (abs(PS_shift_xy(i, 2)) >= 0.1d0*V_max*time_step) then
                    PS_shift_xy(i, 2) = PS_shift_xy(i, 2) / tmp * (0.1d0*V_max*time_step)
                    SP_xy(i, 2) = SP_xy(i, 2) + PS_shift_xy(i, 2)
                else
                    SP_xy(i, 2) = SP_xy(i, 2) + PS_shift_xy(i, 2)
                endif

                ! renew v, rho_S, T
                SP_uv(i, 1) = SP_uv(i, 1) + (PS_TE_uv(i, 1)    * PS_shift_xy(i, 1) + PS_TE_uv(i, 2)    * PS_shift_xy(i, 2))
                SP_uv(i, 2) = SP_uv(i, 2) + (PS_TE_uv(i, 3)    * PS_shift_xy(i, 1) + PS_TE_uv(i, 4)    * PS_shift_xy(i, 2))
                SP_rho_S(i) = SP_rho_S(i) + (PS_TE_rho_S(i, 1) * PS_shift_xy(i, 1) + PS_TE_rho_S(i, 2) * PS_shift_xy(i, 2))
                SP_tem(i)   = SP_tem(i)   + (PS_TE_tem(i, 1)   * PS_shift_xy(i, 1) + PS_TE_tem(i, 2)   * PS_shift_xy(i, 2))
                ! SP_pre(i) = SP_pre(i) + (PS_TE_pre(i, 1) * PS_shift_xy(i, 1) + PS_TE_pre(i, 2) * PS_shift_xy(i, 2))
                ! SP_sm_rho(i) = SP_sm_rho(i) + (PS_TE_sm_rho(i, 1) * PS_shift_xy(i, 1) + PS_TE_sm_rho(i, 2) * PS_shift_xy(i, 2))

                ! rebound
                if (SP_xy(i, 1) < wall_width) then
                    SP_xy(i, 1) = 2.0d0*wall_width - SP_xy(i, 1)
                elseif (SP_xy(i, 1) > (wall_width + length) ) then
                    SP_xy(i, 1) = 2.0d0*(length+wall_width) - SP_xy(i, 1)
                elseif (SP_xy(i, 2) < wall_width) then
                    SP_xy(i, 2) = 2.0d0*wall_width - SP_xy(i, 2)
                elseif (SP_xy(i, 2) > (wall_width + length) ) then
                    SP_xy(i, 2) = 2.0d0*(length+wall_width) - SP_xy(i, 2)
                end if

                ! renew pressure
                call cal_artificial_EOS(SP_rho_S(i), SP_pre(i))

                ! renew rho_G
                call cal_real_EOS(SP_pre(i), SP_tem(i), SP_rho_G(i))
            endif
        enddo
        !$OMP enddo
        !$OMP barrier
        !$OMP end parallel
    end subroutine PS_scheme_exe
end module PS_scheme