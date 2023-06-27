!==========================================================================================!
!                           thermal convection using GDS-WCSPH                             !
!------------------------------------------------------------------------------------------!
!                          Copyright by Kensuke Shobuzako (2022)                           !
!==========================================================================================!

module WRITE_OUT
    
    !$ use omp_lib

    use INPUT
    use GLB_SETTING
    use INITIAL_SETTING
    use analysis

    implicit none

    ! saved
    character(100), save :: tmp_file_name
    real(8), save :: time_loop_start, time_loop_tmp
    integer, save :: write_count
    character(100), save :: save_file_num
    
contains

    !---------------------
    ! program progress
    !---------------------
    subroutine write_progress(t)
        implicit none
        integer, intent(in) :: t
        if ( t == 1 ) then
            tmp_file_name =trim('../output_setting/state_of_progress_'//trim(adjustl(file_name))//'.dat')
            open(50, file=tmp_file_name, status='replace', position='append')
            !$ time_loop_start = omp_get_wtime()
            write(50, '(1x,a)') '--------------------------------------------------------------'
            write(50, '(21x,a)') 'State of progress'
            write(50, '(1x,a)') '--------------------------------------------------------------'
            write(50, '(1x,a,e18.10)') 'Time-loop has been started !!'
            close(50)
        elseif ( t == 10 ) then
            tmp_file_name =trim('../output_setting/state_of_progress_'//trim(adjustl(file_name))//'.dat')
            open(50, file=tmp_file_name, status='old', position='append') ! progress INFO
            !$ time_loop_tmp = omp_get_wtime()
            write(50, '(1x,a,e18.10)') '10 steps have been finished and passed    :' ,(time_loop_tmp)-(time_loop_start)
            close(50)
        elseif ( t == int(total_step/3.0d0) ) then
            tmp_file_name =trim('../output_setting/state_of_progress_'//trim(adjustl(file_name))//'.dat')
            open(50, file=tmp_file_name, status='old', position='append') ! progress INFO
            !$ time_loop_tmp = omp_get_wtime()
            write(50, '(1x,a,e18.10)') '1/3 has been finished and passed          :' ,(time_loop_tmp)-(time_loop_start)
            close(50)
        elseif ( t == 2 * int(total_step/3.0d0) ) then
            tmp_file_name =trim('../output_setting/state_of_progress_'//trim(adjustl(file_name))//'.dat')
            open(50, file=tmp_file_name, status='old', position='append') ! progress INFO
            !$ time_loop_tmp = omp_get_wtime()
            write(50, '(1x,a,e18.10)') '2/3 has been finished and passed          :', (time_loop_tmp)-(time_loop_start)
            close(50)
        elseif ( t == total_step ) then
            tmp_file_name =trim('../output_setting/state_of_progress_'//trim(adjustl(file_name))//'.dat')
            open(50, file=tmp_file_name, status='old', position='append') ! progress INFO
            !$ time_loop_tmp = omp_get_wtime()
            write(50, '(1x,a,e18.10)') 'Time-loop has been completed and passed   :', (time_loop_tmp)-(time_loop_start)
            close(50)
        endif
    end subroutine write_progress


    !---------------------
    ! announce end
    !---------------------
    subroutine write_announce_end(start_time, end_time)
        implicit none
        real(8), intent(in) :: start_time, end_time
        tmp_file_name =trim('../output_setting/state_of_progress_'//trim(adjustl(file_name))//'.dat')
        open(50, file=tmp_file_name, status='old', position='append') ! progress INFO
        write(50, '(1x,a,e18.10)') 'All programs have done, and total time is :', (end_time)-(start_time)
        close(50)

        write(*,*) ''
        write(*,*) '[Message] Your calculation has been completed !!'
        write(*,*) '[Message] calculation time [s] :', (end_time - start_time) 
        write(*,*) ''
    end subroutine write_announce_end


    !---------------------
    ! system setting
    !---------------------
    subroutine write_glb_setting
        ! for system_check
        tmp_file_name =trim('../output_setting/system_check_'//trim(adjustl(file_name))//'.dat')
        open(10, file=tmp_file_name, status='replace')
        write(10, '(1x,a)') '------------------------------------------------'
        write(10, '(13x,a)') 'System INFO (SI unit)'
        write(10, '(1x,a)') '------------------------------------------------'
        write(10, *) 'file_name              :  ', file_name
        write(10, '(1x,a,e18.10)') 'num                    :', dble(num)
        write(10, '(1x,a,e18.10)') 'num_system             :', dble(num_system)
        write(10, '(1x,a,e18.10)') 'total_num              :', dble(total_num)
        write(10, '(1x,a,e18.10)') 'total_num_system       :', dble(total_num_system)
        write(10, '(1x,a,e18.10)') 'total_num_wall         :', dble(total_num_system - total_num)
        write(10, '(1x,a,e18.10)') 'num_MP_side            :', dble(num_MP_side)
        write(10, '(1x,a,e18.10)') 'num_MP                 :', dble(num_MP)
        write(10, '(1x,a,e18.10)') 'wall_thickness         :', dble(wall_thickness)
        write(10, '(1x,a,e18.10)') 'length                 :', dble(length)
        write(10, '(1x,a,e18.10)') 'Delta T                :', delta_tem
        write(10, '(1x,a,e18.10)') 'coe_smooth             :', coe_smooth
        write(10, '(1x,a,e18.10)') 'Delta x                :', delta_x
        write(10, '(1x,a,e18.10)') 'smooth_len             :', smooth_len
        write(10, '(1x,a,e18.10)') 'SP_mass                :', SP_mass
        write(10, '(1x,a,e18.10)') 'sprt_dom               :', sprt_dom
        write(10, '(1x,a,e18.10)') 'zeta                   :', zeta
        write(10, '(1x,a,e18.10)') 'xi                     :', xi
        write(10, '(1x,a,e18.10)') 'coe_step               :', coe_step
        write(10, '(1x,a,e18.10)') 'total_step             :', dble(total_step)
        write(10, '(1x,a,e18.10)') 'writing_step           :', dble(write_step)
        write(10, '(1x,a)') '------------------------------------------------'
        write(10, '(10x,a)') 'Physical Parameters (SI unit)'
        write(10, '(1x,a)') '------------------------------------------------'
        write(10, '(1x,a,e18.10)') 'rho_0                  :', rho_0
        write(10, '(1x,a,e18.10)') 'alpha                  :', alpha
        write(10, '(1x,a,e18.10)') 'c_p                    :', c_p
        write(10, '(1x,a,e18.10)') 'k                      :', k_h
        write(10, '(1x,a,e18.10)') 'eta_0                  :', eta_0
        write(10, '(1x,a,e18.10)') 'K_0                    :', K_0
        write(10, '(1x,a,e18.10)') 'T_0                    :', T_0
        write(10, '(1x,a,e18.10)') 'nu                     :', nu
        write(10, '(1x,a,e18.10)') 'kappa                  :', kappa
        write(10, '(1x,a,e18.10)') 'Ra                     :', Ra
        write(10, '(1x,a,e18.10)') 'Pr                     :', Pr
        write(10, '(1x,a,e18.10)') 'eff_Pr                 :', eff_Pr
        write(10, '(1x,a,e18.10)') 'M_th                   :', M_th
        write(10, '(1x,a,e18.10)') 'SS                     :', SS
        write(10, '(1x,a,e18.10)') 'Reduced SS             :', RSS
        write(10, '(1x,a,e18.10)') 'V_fin                  :', v_fin
        write(10, '(1x,a,e18.10)') 'V_fin_tb               :', v_fin_tb
        write(10, '(1x,a)') '------------------------------------------------'
        write(10, '(15x,a)') 'Time Step (SI unit)'
        write(10, '(1x,a)') '------------------------------------------------'
        write(10, '(1x,a,e18.10)') 'Delta t (CFL)          :', t_CFL
        write(10, '(1x,a,e18.10)') 'Delta t (mo)           :', t_VNM
        write(10, '(1x,a,e18.10)') 'Delta t (th)           :', t_VNT
        write(10, '(1x,a,e18.10)') 'reduced Delta t (CFL)  :', r_t_CFL   
        write(10, '(1x,a,e18.10)') 'reduced Delta t (mo)   :', r_t_VNM
        write(10, '(1x,a,e18.10)') 'time step              :', time_step
        write(10, '(1x,a,e18.10)') 'tau_V_fin              :', tau_V_fin
        write(10, '(1x,a,e18.10)') 'tau_th                 :', tau_th
        write(10, '(1x,a,e18.10)') 'step_V_fin             :', step_V_fin
        write(10, '(1x,a,e18.10)') 'step_th                :', step_th
        !write(10, '(1x,a,e18.10)') 'step_change_s         :', dble(step_change_s)
        !write(10, '(1x,a,e18.10)') 'step_change_e         :', dble(step_change_e)
        write(10, '(1x,a,e18.10)') 'eff_write              :', eff_write
        write(10, '(1x,a,e18.10)') 'total_cal_time         :', total_time
        write(10, '(1x,a)') '------------------------------------------------'
        write(10, '(16x,a)') 'for Python plot'
        write(10, '(1x,a)') '------------------------------------------------'
        write(10, '(a, e16.10, a)') 'total_step    = int(', dble(total_step), ')'
        write(10, '(a, e16.10, a)') 'num           = int(', dble(num), ')'
        write(10, '(a, e16.10, a)') 'num_system    = int(', dble(num_system), ')'
        write(10, '(a, e16.10, a)') 'total_num     = int(', dble(total_num), ')'
        write(10, '(a, e16.10, a)') 'total_num_system = int(', dble(total_num_system), ')'
        write(10, '(a, e16.10, a)') 'num_wall      = int(', dble(total_num_system - total_num), ')'
        write(10, '(a, e16.10, a)') 'num_MP_side   = int(', dble(num_MP_side), ')'
        write(10, '(a, e16.10, a)') 'num_MP        = int(', dble(num_MP), ')'
        write(10, '(a, e16.10, a)') 'wall_thickness = int(', dble(wall_thickness), ')'
        write(10, '(a, e16.10, a)') 'write_step    = int(', dble(write_step), ')'
        write(10, '(a, e16.10, a)') 'RSST_VIM      = int(', dble(RSST_VIM), ')'
        write(10, '(a, e16.10, a)') 'kernel_switch = int(', dble(kernel_switch), ')'
        write(10, *) ''
        write(10, '(a, e16.10)') 'length         = ', length
        write(10, '(a, e16.10)') 'delta_tem      = ', delta_tem
        write(10, '(a, e16.10)') 'coe_smooth     = ', coe_smooth
        write(10, '(a, e16.10)') 'delta_x        = ', delta_x
        write(10, '(a, e16.10)') 'smooth_len     = ', smooth_len
        write(10, '(a, e16.10)') 'support_domain = ', sprt_dom
        write(10, '(a, e16.10)') 'time_step      = ', time_step
        write(10, *) ''
        write(10, '(a, e16.10)') 'rho_0   = ', rho_0
        write(10, '(a, e16.10)') 'alpha   = ', alpha
        write(10, '(a, e16.10)') 'eta_0   = ', eta_0
        write(10, '(a, e16.10)') 'T_0     = ', T_0
        write(10, '(a, e16.10)') 'T_top   = ', T_top
        write(10, '(a, e16.10)') 'kappa   = ', kappa
        write(10, *) ''
        write(10, '(a, e16.10)') 'V_fin    = ', V_fin
        write(10, '(a, e16.10)') 'V_fin_tb = ', V_fin_tb
        write(10, *) ''
        write(10, '(a, e16.10)') 'Ra   = ', Ra
        write(10, '(a, e16.10)') 'Pr   = ', Pr
        write(10, '(a, e16.10)') 'M_th = ', M_th
        write(10, '(a, e16.10)') 'SS   = ', SS
        write(10, '(a, e16.10)') 'RSS  = ', RSS
        close(10)

        ! for python
        tmp_file_name =trim('../output_setting/for_python_'//trim(adjustl(file_name))//'.dat')
        open(11, file=tmp_file_name, status='replace')
        write(11, '(e16.10)') dble(total_step)
        write(11, '(e16.10)') dble(num)
        write(11, '(e16.10)') dble(num_system)
        write(11, '(e16.10)') dble(total_num)
        write(11, '(e16.10)') dble(total_num_system)
        write(11, '(e16.10)') dble(total_num_system - total_num)
        write(11, '(e16.10)') dble(num_MP_side)
        write(11, '(e16.10)') dble(num_MP)
        write(11, '(e16.10)') dble(wall_thickness)
        write(11, '(e16.10)') dble(write_step)
        write(11, '(e16.10)') dble(RSST_VIM)
        write(11, '(e16.10)') dble(kernel_switch)
        ! write(11, *) ''
        write(11, '(e16.10)') length
        write(11, '(e16.10)') delta_tem
        write(11, '(e16.10)') coe_smooth
        write(11, '(e16.10)') delta_x
        write(11, '(e16.10)') smooth_len
        write(11, '(e16.10)') sprt_dom
        write(11, '(e16.10)') time_step
        ! write(11, *) ''
        write(11, '(e16.10)') rho_0
        write(11, '(e16.10)') alpha
        write(11, '(e16.10)') eta_0
        write(11, '(e16.10)') T_0
        write(11, '(e16.10)') T_top
        write(11, '(e16.10)') kappa
        ! write(11, *) ''
        write(11, '(e16.10)') V_fin
        write(11, '(e16.10)') V_fin_tb
        ! write(11, *) ''
        write(11, '(e16.10)') Ra
        write(11, '(e16.10)') Pr
        write(11, '(e16.10)') M_th
        write(11, '(e16.10)') SS
        write(11, '(e16.10)') RSS
        close(11)

    end subroutine write_glb_setting


    !---------------------
    ! initial setting
    !---------------------
    subroutine write_initial_setting
        ! initial state
        tmp_file_name = trim('../output_setting/initial_state_SP_'//trim(adjustl(file_name))//'.dat')
        open(20, file=tmp_file_name, status='replace', form='unformatted', access='stream')
        write(20) SP_xy, SP_uv, SP_rho_S, SP_pre, SP_tem, SP_eta, SP_sm_rho, SP_rho_G
        close(20)

        ! particle ID
        tmp_file_name = trim('../output_setting/initial_state_SP_kind_'//trim(adjustl(file_name))//'.dat')
        open(21, file=tmp_file_name, status='replace', form='unformatted', access='stream')
        write(21) SP_kind
        close(21)

        ! NNPS
        tmp_file_name = trim('../output_setting/initial_state_NNPS_'//trim(adjustl(file_name))//'.dat')
        open(22, file=tmp_file_name, status='replace', form='unformatted', access='stream')
        write(22) int(SP_uv(:, 1))
        close(22)
        SP_uv(:, :) = 0.0d0

        ! VM
        tmp_file_name = trim('../output_setting/initial_VM_'//trim(adjustl(file_name))//'.dat')
        open(23, file=tmp_file_name, status='replace', form='unformatted', access='stream')
        write(23) VM_xy, VM_uv, VM_pre, VM_tem, ID_VM_WL
        close(23)

        ! NNPS for VM
        tmp_file_name = trim('../output_setting/NNPS_VM_'//trim(adjustl(file_name))//'.dat')
        open(24, file=tmp_file_name, status='replace', form='unformatted', access='stream')
        write(24) int(VM_uv(:, 1))
        close(24)
        VM_uv(:, :) = 0.0d0

        ! MP
        tmp_file_name = trim('../output_setting/initial_MP_'//trim(adjustl(file_name))//'.dat')
        open(25, file=tmp_file_name, status='replace', form='unformatted', access='stream')
        write(25) MP_xy, MP_uv, MP_pre, MP_tem, MP_NU
        close(25)

        ! NNPS for MP
        tmp_file_name = trim('../output_setting/NNPS_MP_'//trim(adjustl(file_name))//'.dat')
        open(27, file=tmp_file_name, status='replace', form='unformatted', access='stream')
        write(27) int(MP_uv(:, 1))
        close(27)
        MP_uv(:, :) = 0.0d0
    end subroutine write_initial_setting

    !---------------------
    ! make at each time
    !---------------------
    subroutine write_save_file
        write(save_file_num, *) write_count
        ! SP
        tmp_file_name = trim('../save_file/SP_save/SP_save_'//trim(adjustl(file_name))// &
                                                         '_'//trim(adjustl(save_file_num))//'.dat')
        open(20, file=tmp_file_name, status='replace', form='unformatted', access='stream')
        write(20) SP_xy, SP_uv, SP_rho_S, SP_pre, SP_tem, SP_eta, SP_sm_rho, SP_rho_G, par_con
        close(20)

        ! EOM
        tmp_file_name = trim('../save_file/EOM_save/EOM_save_'//trim(adjustl(file_name))// &
                                                         '_'//trim(adjustl(save_file_num))//'.dat')
        open(21, file=tmp_file_name, status='replace', form='unformatted', access='stream')
        write(21) vis_EOM, buo_EOM, pre_EOM, EOM_b*xi**(2.0d0)
        close(21)
            
        ! MP
        tmp_file_name = trim('../save_file/MP_save/MP_save_'//trim(adjustl(file_name))//&
                                                 '_'//trim(adjustl(save_file_num))//'.dat')
        open(30, file=tmp_file_name, status='replace', form='unformatted', access='stream')
        write(30) MP_xy, MP_uv, MP_pre, MP_tem, MP_NU
        close(30)
            
        ! V_rms
        if (t == 0) then
            tmp_file_name = trim('../save_file/V_rms_save_'//trim(adjustl(file_name))//'.dat')
            open(40, file=tmp_file_name, status='replace')
            write(40, '(e18.10)') V_rms_ave
            close(40)
        else
            tmp_file_name = trim('../save_file/V_rms_save_'//trim(adjustl(file_name))//'.dat')
            open(40, file=tmp_file_name, status='old', position='append')
            write(40, '(e18.10)') V_rms_ave
            close(40)
        endif

        ! time_save
        if (t == 0) then
            tmp_file_name = trim('../save_file/time_save_'//trim(adjustl(file_name))//'.dat')
            open(41, file=tmp_file_name, status='replace')
            write(41, '(e18.10)') pass_time
            close(41)
        else
            tmp_file_name = trim('../save_file/time_save_'//trim(adjustl(file_name))//'.dat')
            open(41, file=tmp_file_name, status='old', position='append')
            write(41, '(e18.10)') pass_time
            close(41)
        endif

        ! ! delta_t
        ! if (t == 0) then
        !     tmp_file_name = trim('../save_file/time_step_'//trim(adjustl(file_name))//'.dat')
        !     open(42, file=tmp_file_name, status='replace')
        !     write(42, '(e18.10)') time_step
        !     close(42)
        ! else
        !     tmp_file_name = trim('../save_file/time_step_'//trim(adjustl(file_name))//'.dat')
        !     open(42, file=tmp_file_name, status='old', position='append')
        !     write(42, '(e18.10)') time_step
        !     close(42)
        ! endif

        ! ! RSS
        ! if (t == 0) then
        !     tmp_file_name = trim('../save_file/RSS_'//trim(adjustl(file_name))//'.dat')
        !     open(43, file=tmp_file_name, status='replace')
        !     write(43, '(e18.10)') RSS
        !     close(43)
        ! else
        !     tmp_file_name = trim('../save_file/RSS_'//trim(adjustl(file_name))//'.dat')
        !     open(43, file=tmp_file_name, status='old', position='append')
        !     write(43, '(e18.10)') RSS
        !     close(43)
        ! endif

        ! PS_scheme
        tmp_file_name = trim('../save_file/PS_scheme_save/PS_scheme_save_'&
                            //trim(adjustl(file_name))//'_'//trim(adjustl(save_file_num))//'.dat')
        open(50, file=tmp_file_name, status='replace', form='unformatted', access='stream')
        write(50) PS_shift_xy
        close(50)

        write_count = write_count + 1
    end subroutine write_save_file

end module WRITE_OUT