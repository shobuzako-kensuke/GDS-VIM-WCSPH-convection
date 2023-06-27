!==========================================================================================!
!                           thermal convection using GDS-WCSPH                             !
!------------------------------------------------------------------------------------------!
!                          Copyright by Kensuke Shobuzako (2022)                           !
!==========================================================================================!

module cal_WALL

    use INPUT
    use GLB_SETTING
    use GLB_ARG
    use INITIAL_SETTING
    use NNPS
    use cal_KERNEL

    implicit none
    integer :: my_cell, your_cell
    real(8) :: x_ij, y_ij, harm_eta, vis_x, vis_y, u_ij, v_ij
    real(8) :: tmp_W_ij, tmp_thick

contains
    subroutine cal_VM_WALL
        !-----------------------
        ! F_i_and_sm_rho
        !-----------------------
        ! zero clear
        F_i(:,:)        = 0.0d0
        SP_sm_rho(:)    = 0.0d0
        par_con(:)      = 0.0d0
        VM_corrected(:) = 0.0d0
        VM_uv(:,:)      = 0.0d0
        VM_pre(:)       = 0.0d0
        VM_tem(:)       = 0.0d0
        VM_rho(:)       = 0.0d0
        VM_F_i(:,:)     = 0.0d0
        
        !$ call omp_set_dynamic(.false.)
        !$ call omp_set_num_threads(num_threads)
        !$OMP parallel default(none) &
        !$OMP private(i, j, k, me, you, my_cell, your_cell) &
        !$OMP private(x_ij, y_ij, r_ij, q_ij, W_ij, dW_ij) &
        !$OMP private(tmp, u_ij, v_ij, harm_eta, vis_x, vis_y, tmp_W_ij) &
        !$OMP shared(inner_cell_start, inner_cell_end, total_num_system) &
        !$OMP shared(NNPS_arg, SP_kind, NNPS_9, NNPS_max, smooth_len) &
        !$OMP shared(SP_mass, SP_xy, SP_uv, SP_rho_S, SP_pre) &
        !$OMP shared(SP_tem, SP_eta, SP_rho_G, SP_sm_rho, F_i) &
        !$OMP shared(NNPS_row, NNPS_VM_row) &
        !$OMP shared(VM_corrected, NNPS_arg_VM) &
        !$OMP shared(VM_xy, VM_uv, VM_pre, VM_tem, VM_rho, VM_F_i) &
        !$OMP shared(num_VM, ID_VM_WL, wall_thickness, delta_x, rho_0, length) &
        !$OMP shared(T_bot, T_top, par_con, wall_width)

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

                            !-----------------------
                            ! cal kernel
                            !-----------------------
                            x_ij = SP_xy(me, 1) - SP_xy(you, 1)
                            y_ij = SP_xy(me, 2) - SP_xy(you, 2)
                            r_ij = dsqrt(x_ij*x_ij + y_ij*y_ij)
                            q_ij = r_ij / smooth_len
                            call cal_W_ij_and_dW_ij(q_ij, W_ij, dW_ij)

                            !-----------------------
                            ! cal EOM
                            !-----------------------
                            if (me /= you) then
                                tmp = SP_mass / (SP_rho_S(me) * SP_rho_S(you))
                                u_ij = SP_uv(me, 1) - SP_uv(you, 1)
                                v_ij = SP_uv(me, 2) - SP_uv(you, 2)

                                harm_eta = 4.0d0 * SP_eta(me) * SP_eta(you) / (SP_eta(me) + SP_eta(you))
                                vis_x = harm_eta * tmp * dW_ij / r_ij * u_ij
                                vis_y = harm_eta * tmp * dW_ij / r_ij * v_ij

                                F_i(me, 1) = F_i(me, 1) + vis_x
                                F_i(me, 2) = F_i(me, 2) + vis_y
                            endif

                            !-----------------------
                            ! cal smooth_rho
                            !-----------------------
                            SP_sm_rho(me) = SP_sm_rho(me) + SP_mass * W_ij

                            !-----------------------
                            ! cal par_con
                            !-----------------------
                            par_con(me) = par_con(me) + (SP_mass / SP_rho_S(you)) * W_ij
                        enddo
                    enddo
                endif
            enddo
        enddo
        !$OMP enddo
        !$OMP barrier

        ! F_i(:,2) += gravity
        !$OMP do
        do i = 1, total_num_system
            F_i(i, 2) = F_i(i, 2) + SP_rho_G(i) / SP_rho_S(i) * g
        enddo
        !$OMP enddo
        !$OMP barrier

        !-----------------------
        ! corrected W for VM
        !-----------------------
        !$OMP do
        ! my cell number loop
        do my_cell = inner_cell_start, inner_cell_end

            ! my number loop
            do i = 1, NNPS_VM_row(my_cell)-1

                ! if no VM exists, exit
                if (NNPS_VM_row(my_cell) == 1) exit

                ! my number
                me = NNPS_arg_VM(i, my_cell)

                ! your cell loop
                do j = 1, 9
                    your_cell = NNPS_9(j, my_cell)

                    ! your number loop
                    do k = 1, NNPS_row(your_cell)-1
                        you = NNPS_arg(k, your_cell)

                        ! if you are inner_particle, exe
                        if (SP_kind(you) == 0) then
                            !-----------------------
                            ! cal kernel
                            !-----------------------
                            x_ij = VM_xy(me, 1) - SP_xy(you, 1)
                            y_ij = VM_xy(me, 2) - SP_xy(you, 2)
                            r_ij = dsqrt(x_ij*x_ij + y_ij*y_ij)
                            q_ij = r_ij / smooth_len
                            call cal_W_ij(q_ij, W_ij)

                            !-----------------------
                            ! cal corrected W
                            !-----------------------
                            tmp = SP_mass / SP_rho_S(you)
                            VM_corrected(me) = VM_corrected(me) + tmp * W_ij
                        endif
                    enddo
                enddo
            enddo
        enddo
        !$OMP enddo
        !$OMP barrier

        !-----------------------
        ! cal_VM
        !-----------------------
        !$OMP do
        ! my cell number loop
        do my_cell = inner_cell_start, inner_cell_end

            ! my number loop
            do i = 1, NNPS_VM_row(my_cell)-1

                ! if no VM exists, exit
                if (NNPS_VM_row(my_cell) == 1) exit

                ! my number
                me = NNPS_arg_VM(i, my_cell)

                ! your cell loop
                do j = 1, 9
                    your_cell = NNPS_9(j, my_cell)

                    ! your number loop
                    do k = 1, NNPS_row(your_cell)-1
                        you = NNPS_arg(k, your_cell)

                        ! if you are inner_particle, exe
                        if (SP_kind(you) == 0) then
                            !-----------------------
                            ! cal kernel
                            !-----------------------
                            x_ij = VM_xy(me, 1) - SP_xy(you, 1)
                            y_ij = VM_xy(me, 2) - SP_xy(you, 2)
                            r_ij = dsqrt(x_ij*x_ij + y_ij*y_ij)
                            q_ij = r_ij / smooth_len
                            call cal_W_ij(q_ij, W_ij)

                            !-----------------------
                            ! cal correted VM
                            !-----------------------
                            tmp = SP_mass / SP_rho_S(you)
                            if (VM_corrected(me) == 0.0d0) then
                                VM_corrected(me) = 1.0d0
                            endif
                            tmp_W_ij = W_ij / VM_corrected(me)

                            VM_uv(me, 1)  = VM_uv(me, 1)  + SP_uv(you, 1) * tmp * tmp_W_ij
                            VM_uv(me, 2)  = VM_uv(me, 2)  + SP_uv(you, 2) * tmp * tmp_W_ij
                            VM_pre(me)    = VM_pre(me)    + SP_pre(you)   * tmp * tmp_W_ij
                            VM_tem(me)    = VM_tem(me)    + SP_tem(you)   * tmp * tmp_W_ij
                            VM_rho(me)    = VM_rho(me)    + SP_rho_S(you) * tmp * tmp_W_ij
                            VM_F_i(me, 1) = VM_F_i(me, 1) + F_i(you, 1)   * tmp * tmp_W_ij
                            VM_F_i(me, 2) = VM_F_i(me, 2) + F_i(you, 2)   * tmp * tmp_W_ij
                        endif
                    enddo
                enddo
            enddo
        enddo
        !$OMP enddo
        !$OMP barrier

        !-----------------------
        ! VM -> wall
        !-----------------------
        !$OMP do
        do me = 1, num_VM
            you = ID_VM_WL(me)
            x_ij = VM_xy(me, 1) - SP_xy(you, 1)
            y_ij = VM_xy(me, 2) - SP_xy(you, 2)

            ! bottom wall
            if (SP_kind(you) == 1) then
                ! bottom right or bottom left wall
                if ((SP_xy(you, 1) < wall_width) .or. (SP_xy(you, 1) > (wall_width + length))) then
                    ! slip condition
                    SP_uv(you, 1) = - VM_uv(me, 1)
                    SP_uv(you, 2) = - VM_uv(me, 2)
                ! bottom wall
                else
                    ! slip condition
                    if (wall_v == 0) then
                        SP_uv(you, 1) = - VM_uv(me, 1)
                        SP_uv(you, 2) = - VM_uv(me, 2)
                    else
                        SP_uv(you, 1) = + VM_uv(me, 1)
                        SP_uv(you, 2) = - VM_uv(me, 2)
                    endif
                endif
                SP_pre(you) = VM_pre(me) - y_ij * VM_rho(me) * VM_F_i(me, 2)
                SP_tem(you) = T_bot
            
            ! top wall
            elseif (SP_kind(you) == 2) then
                ! top right or top left wall
                if ((SP_xy(you, 1) < wall_width) .or. (SP_xy(you, 1) > (wall_width + length))) then
                    ! slip condition
                    SP_uv(you, 1) = - VM_uv(me, 1)
                    SP_uv(you, 2) = - VM_uv(me, 2)
                ! top wall
                else
                    !slip condition
                    if (wall_v == 0) then
                        SP_uv(you, 1) = - VM_uv(me, 1)
                        SP_uv(you, 2) = - VM_uv(me, 2)
                    else
                        SP_uv(you, 1) = + VM_uv(me, 1)
                        SP_uv(you, 2) = - VM_uv(me, 2)
                    endif
                endif
                SP_pre(you) = VM_pre(me) - y_ij * VM_rho(me) * VM_F_i(me, 2)
                SP_tem(you) = T_top

            ! right & left wall
            else
                ! slip condition
                if (wall_v == 0) then
                    SP_uv(you, 1) = - VM_uv(me, 1)
                    SP_uv(you, 2) = - VM_uv(me, 2)
                else
                    SP_uv(you, 1) = - VM_uv(me, 1)
                    SP_uv(you, 2) = + VM_uv(me, 2)
                endif
                SP_pre(you) = VM_pre(me) - x_ij * VM_rho(me) * VM_F_i(me, 1)
                SP_tem(you) = VM_tem(me)
            endif
        enddo
        !---------------------------------------------
        ! メモ：粘性率変える場合は壁も変えた方が良いかも
        !---------------------------------------------
        !$OMP enddo
        !$OMP end parallel
    end subroutine cal_VM_WALL
end module cal_WALL