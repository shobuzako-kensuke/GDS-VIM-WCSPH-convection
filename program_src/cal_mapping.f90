!==========================================================================================!
!                           thermal convection using GDS-WCSPH                             !
!------------------------------------------------------------------------------------------!
!                          Copyright by Kensuke Shobuzako (2022)                           !
!==========================================================================================!

module cal_MAPPING
    use INPUT
    use GLB_SETTING
    use GLB_ARG
    use INITIAL_SETTING
    use NNPS
    use cal_KERNEL

    implicit none
    integer :: my_cell, your_cell
    real(8) :: x_ij, y_ij, dW_ij_x, dW_ij_y, tmp_W_ij, T_ji

contains
    subroutine interpolate_MP
        ! zero clear
        MP_corrected(:) = 0.0d0
        MP_KGC(:,:)     = 0.0d0
        MP_uv(:,:)      = 0.0d0
        MP_pre(:)       = 0.0d0
        MP_tem(:)       = 0.0d0
        MP_NU(:)        = 0.0d0

        !$ call omp_set_dynamic(.false.)
        !$ call omp_set_num_threads(num_threads)
        !$OMP parallel default(none) &
        !$OMP private(i, j, k, me, you, my_cell, your_cell) &
        !$OMP private(x_ij, y_ij, r_ij, q_ij, W_ij, dW_ij) &
        !$OMP private(dW_ij_x, dW_ij_y, tmp, tmp_W_ij, T_ji) &
        !$OMP shared(NNPS_max, NNPS_arg_MP, NNPS_9, NNPS_arg) &
        !$OMP shared(SP_xy, SP_uv, SP_rho_S, SP_pre, SP_tem) &
        !$OMP shared(smooth_len, SP_mass, MP_corrected, MP_KGC) &
        !$OMP shared(MP_xy, MP_uv, MP_pre, MP_tem, MP_NU) &
        !$OMP shared(wall_thickness, delta_x, T_bot, T_top, delta_tem, length) &
        !$OMP shared(SP_kind, NNPS_row, wall_width)

        !-----------------------
        ! cal right-hands
        !-----------------------
        !$OMP do
        ! my cell number loop
        do my_cell = 1, NNPS_max

            ! my number loop
            do i = 1, size(NNPS_arg_MP(:, 1))
                me = NNPS_arg_MP(i, my_cell) ! my number

                ! if no MP exists, exit
                if (me == 0) exit

                ! your cell loop
                do j = 1, 9
                    your_cell = NNPS_9(j, my_cell)

                    ! if out of range, exit
                    if ((your_cell <= 0) .or. (NNPS_max < your_cell)) exit 
                    
                    ! your number loop
                    do k = 1, NNPS_row(your_cell)-1
                        you = NNPS_arg(k, your_cell)

                        !-----------------------
                        ! cal kernel
                        !-----------------------
                        x_ij = MP_xy(me, 1) - SP_xy(you, 1)
                        y_ij = MP_xy(me, 2) - SP_xy(you, 2)
                        r_ij = dsqrt(x_ij*x_ij + y_ij*y_ij)
                        q_ij = r_ij / smooth_len
                        call cal_W_ij_and_dW_ij(q_ij, W_ij, dW_ij)

                        !-----------------------
                        ! cal KGC arg
                        !-----------------------
                        tmp = SP_mass / SP_rho_S(you)
                        MP_corrected(me) = MP_corrected(me) + tmp * W_ij

                        ! if you are inner_particle, exe
                        if (SP_kind(you) == 0) then
                            dW_ij_x = x_ij / r_ij * dW_ij
                            dW_ij_y = y_ij / r_ij * dW_ij
                            MP_KGC(me, 1) = MP_KGC(me, 1) + tmp * (-x_ij) * dW_ij_x
                            MP_KGC(me, 2) = MP_KGC(me, 2) + tmp * (-y_ij) * dW_ij_x
                            MP_KGC(me, 3) = MP_KGC(me, 3) + tmp * (-x_ij) * dW_ij_y
                            MP_KGC(me, 4) = MP_KGC(me, 4) + tmp * (-y_ij) * dW_ij_y
                        endif
                    enddo
                enddo
            enddo
        enddo
        !$OMP enddo
        !$OMP barrier

        !-----------------------
        ! cal inverse KGC arg
        !-----------------------
        !$OMP do
        do i = 1, num_MP
            tmp = 1.0d0 / (MP_KGC(i, 1)*MP_KGC(i, 4) - MP_KGC(i, 2)*MP_KGC(i, 3))
            x_ij = MP_KGC(i, 1)
            ! MP_KGC(i, 1) =   tmp * MP_KGC(i, 4)
            ! MP_KGC(i, 2) = - tmp * MP_KGC(i, 2)
            MP_KGC(i, 3) = - tmp * MP_KGC(i, 3)
            MP_KGC(i, 4) =   tmp * x_ij
        enddo
        !$OMP enddo
        !$OMP barrier

        !-----------------------
        ! MP interpolating
        !-----------------------
        !$OMP do
        ! my cell number loop
        do my_cell = 1, NNPS_max 

            ! my number loop
            do i = 1, size(NNPS_arg_MP(:, 1))
                me = NNPS_arg_MP(i, my_cell) ! my number

                ! if not exist NNPS_arg, exit
                if (me == 0) exit

                ! your cell loop
                do j = 1, 9
                    your_cell = NNPS_9(j, my_cell)

                    ! if out of range, exit
                    if ((your_cell <= 0) .or. (NNPS_max < your_cell)) exit 
                    
                    ! your number loop
                    do k = 1, NNPS_row(your_cell)-1
                        you = NNPS_arg(k, your_cell)

                        !-----------------------
                        ! cal kernel
                        !-----------------------
                        x_ij = MP_xy(me, 1) - SP_xy(you, 1)
                        y_ij = MP_xy(me, 2) - SP_xy(you, 2)
                        r_ij = dsqrt(x_ij*x_ij + y_ij*y_ij)
                        q_ij = r_ij / smooth_len
                        call cal_W_ij_and_dW_ij(q_ij, W_ij, dW_ij)

                        !-----------------------
                        ! interpolate
                        !-----------------------
                        tmp = SP_mass / SP_rho_S(you)
                        if (MP_corrected(me) == 0.0d0) then
                            MP_corrected(me) = 1.0d0
                        endif
                        tmp_W_ij = W_ij / MP_corrected(me)

                        MP_uv(me, 1) = MP_uv(me, 1) + SP_uv(you, 1) * tmp * tmp_W_ij
                        MP_uv(me, 2) = MP_uv(me, 2) + SP_uv(you, 2) * tmp * tmp_W_ij
                        MP_pre(me)   = MP_pre(me)   + SP_pre(you)   * tmp * tmp_W_ij
                        MP_tem(me)   = MP_tem(me)   + SP_tem(you)   * tmp * tmp_W_ij
                    enddo
                enddo
            enddo
        enddo
        !$OMP enddo

        !-----------------------
        ! cal Nu number
        !-----------------------
        !$OMP do
        ! my cell number loop
        do my_cell = 1, NNPS_max 

            ! my number loop
            do i = 1, size(NNPS_arg_MP(:, 1))
                me = NNPS_arg_MP(i, my_cell) ! my number

                ! if not exist NNPS_arg, exit
                if (me == 0) exit

                ! your cell loop
                do j = 1, 9
                    your_cell = NNPS_9(j, my_cell)

                    ! if out of range, exit
                    if ((your_cell <= 0) .or. (NNPS_max < your_cell)) exit 
                    
                    ! your number loop
                    do k = 1, NNPS_row(your_cell)-1
                        you = NNPS_arg(k, your_cell)

                        ! if you are inner_particle, exe
                        if (SP_kind(you) == 0) then
                            !-----------------------
                            ! cal kernel
                            !-----------------------
                            x_ij = MP_xy(me, 1) - SP_xy(you, 1)
                            y_ij = MP_xy(me, 2) - SP_xy(you, 2)
                            r_ij = dsqrt(x_ij*x_ij + y_ij*y_ij)
                            q_ij = r_ij / smooth_len
                            call cal_dW_ij(q_ij, dW_ij)

                            !-----------------------
                            ! interpolate
                            !-----------------------
                            tmp = SP_mass / SP_rho_S(you)
                            dW_ij_x = x_ij / r_ij * dW_ij
                            dW_ij_y = y_ij / r_ij * dW_ij
                            T_ji = SP_tem(you) - MP_tem(me)
                            MP_NU(me) = MP_NU(me) + tmp * T_ji * (MP_KGC(me, 3)*dW_ij_x + MP_KGC(me, 4)*dW_ij_y) &
                                                    / (delta_tem / length)
                            !MP_NU(me) = MP_NU(me) + tmp * T_ji * dW_ij_y / (delta_tem / length)
                        endif
                    enddo
                enddo
            enddo
        enddo
        !$OMP enddo
        !$OMP end parallel
    end subroutine interpolate_MP

end module cal_MAPPING