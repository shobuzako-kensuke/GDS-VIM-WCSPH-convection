!==========================================================================================!
!                           thermal convection using GDS-WCSPH                             !
!------------------------------------------------------------------------------------------!
!                          Copyright by Kensuke Shobuzako (2022)                           !
!==========================================================================================!

module cal_MAIN

    use INPUT
    use GLB_SETTING
    use cal_EOS
    use INITIAL_SETTING
    use NNPS
    use cal_KERNEL
    use set_time_step

    implicit none
    integer :: my_cell, your_cell
    real(8) :: x_ij, y_ij, dW_ij_x, dW_ij_y, u_ij, v_ij
    real(8) :: pre_x, pre_y, vis_x, vis_y, harm_eta
    real(8) :: T_ij

contains

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! 2nd-order Adams-Bashforth
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine cal_main_2AB
        ! zero clear
        vis_EOM(:,:) = 0.0d0
        buo_EOM(:)   = 0.0d0
        pre_EOM(:,:) = 0.0d0

        !$ call omp_set_dynamic(.false.)
        !$ call omp_set_num_threads(num_threads)
        !$OMP parallel default(none) &
        !$OMP private(i, j, k, me, you, my_cell, your_cell) &
        !$OMP private(x_ij, y_ij, r_ij, q_ij, dW_ij) &
        !$OMP private(dW_ij_x, dW_ij_y, u_ij, v_ij) &
        !$OMP private(tmp, harm_eta, pre_x, pre_y, vis_x, vis_y) &
        !$OMP private(T_ij) &
        !$OMP shared(inner_cell_start, inner_cell_end) &
        !$OMP shared(NNPS_max, total_num_system, NNPS_arg, SP_kind, NNPS_9) &
        !$OMP shared(SP_mass, SP_xy, SP_uv, SP_rho_S, SP_pre, SP_tem, SP_eta, SP_rho_G) &
        !$OMP shared(EOC_b, EOM_b, EOE_b, EOC_n, EOM_n, EOE_n) &
        !$OMP shared(k_h, c_p, smooth_len, length, delta_x, wall_thickness) &
        !$OMP shared(t, time_step, wall_width, NNPS_row) &
        !$OMP shared(vis_EOM, buo_EOM, pre_EOM)

        !-----------------------
        ! cal right-hands
        !-----------------------
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
                                ! cal EOC
                                !-----------------------
                                dW_ij_x = x_ij / r_ij * dW_ij
                                dW_ij_y = y_ij / r_ij * dW_ij
                                u_ij = SP_uv(me, 1) - SP_uv(you, 1)
                                v_ij = SP_uv(me, 2) - SP_uv(you, 2)
                                    
                                EOC_n(me) = EOC_n(me) + SP_mass * (u_ij * dW_ij_x + v_ij * dW_ij_y)

                                !-----------------------
                                ! cal EOM
                                !-----------------------
                                tmp = SP_mass / (SP_rho_S(me) * SP_rho_S(you))
                                harm_eta = 4.0d0 * SP_eta(me) * SP_eta(you) / (SP_eta(me) + SP_eta(you))

                                pre_EOM(me, 1) = pre_EOM(me, 1) - tmp * (SP_pre(me) + SP_pre(you)) * dW_ij_x
                                pre_EOM(me, 2) = pre_EOM(me, 2) - tmp * (SP_pre(me) + SP_pre(you)) * dW_ij_y
                                vis_EOM(me, 1) = vis_EOM(me, 1) + harm_eta * tmp * dW_ij / r_ij * u_ij
                                vis_EOM(me, 2) = vis_EOM(me, 2) + harm_eta * tmp * dW_ij / r_ij * v_ij

                                !-----------------------
                                ! cal EOE
                                !-----------------------
                                T_ij = SP_tem(me) - SP_tem(you)

                                EOE_n(me) = EOE_n(me) + (2.0d0*k_h / c_p) * tmp * dW_ij / r_ij * T_ij
                            endif
                        enddo
                    enddo
                endif
            enddo
        enddo
        !$OMP enddo
        !$OMP barrier

        !-----------------------
        ! renewal
        !-----------------------
        !$OMP do
        do i = 1, total_num_system
            if (SP_kind(i) == 0) then
                ! right hand of EOM
                buo_EOM(i) = SP_rho_G(i) / SP_rho_S(i) * g
                EOM_n(i, 1) = pre_EOM(i, 1) + vis_EOM(i, 1)
                EOM_n(i, 2) = pre_EOM(i, 2) + vis_EOM(i, 2) + buo_EOM(i)

                ! multiply xi
                EOM_n(i, :) = EOM_n(i ,:) / (xi*xi)

                ! renew rho_S, v, T, x
                if (t == 1) then
                    SP_rho_S(i) = SP_rho_S(i) + time_step * EOC_n(i)
                    SP_uv(i, 1) = SP_uv(i, 1) + time_step * EOM_n(i, 1)
                    SP_uv(i, 2) = SP_uv(i, 2) + time_step * EOM_n(i, 2)
                    SP_tem(i)   = SP_tem(i)   + time_step * EOE_n(i)
                    SP_xy(i, 1) = SP_xy(i, 1) + time_step * SP_uv(i, 1)
                    SP_xy(i, 2) = SP_xy(i, 2) + time_step * SP_uv(i, 2)
                else
                    SP_rho_S(i) = SP_rho_S(i) + (1.5d0 * time_step * EOC_n(i)    - 0.5d0 * time_step * EOC_b(i))
                    SP_uv(i, 1) = SP_uv(i, 1) + (1.5d0 * time_step * EOM_n(i, 1) - 0.5d0 * time_step * EOM_b(i, 1))
                    SP_uv(i, 2) = SP_uv(i, 2) + (1.5d0 * time_step * EOM_n(i, 2) - 0.5d0 * time_step * EOM_b(i, 2))
                    SP_tem(i)   = SP_tem(i)   + (1.5d0 * time_step * EOE_n(i)    - 0.5d0 * time_step * EOE_b(i))
                    SP_xy(i, 1) = SP_xy(i, 1) + time_step * SP_uv(i, 1)
                    SP_xy(i, 2) = SP_xy(i, 2) + time_step * SP_uv(i, 2)
                endif

                ! renew pressure
                call cal_artificial_EOS(SP_rho_S(i), SP_pre(i))

                ! renew rho_G
                call cal_real_EOS(SP_pre(i), SP_tem(i), SP_rho_G(i))

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
            endif
        enddo
        !$OMP enddo
        !$OMP barrier

        !------------------------------
        ! メモ：粘性率変化を後で追加する
        !------------------------------

        !$OMP end parallel

        ! renewal and clear
        EOC_b(:)   = EOC_n(:)
        EOM_b(:,:) = EOM_n(:,:)
        EOE_b(:)   = EOE_n(:)
        EOC_n(:)   = 0.0d0
        EOM_n(:,:) = 0.0d0
        EOE_n(:)   = 0.0d0
    end subroutine cal_main_2AB
    
end module cal_MAIN