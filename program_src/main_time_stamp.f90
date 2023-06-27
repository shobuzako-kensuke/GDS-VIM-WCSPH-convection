!==========================================================================================!
!                           thermal convection using GDS-WCSPH                             !
!------------------------------------------------------------------------------------------!
!                          Copyright by Kensuke Shobuzako (2022)                           !
!==========================================================================================!

program main
    use INPUT
    use GLB_SETTING
    use GLB_ARG
    use cal_EOS
    use INITIAL_SETTING
    use NNPS
    use cal_KERNEL
    use cal_WALL
    use cal_MAIN
    use PS_scheme
    use cal_MAPPING
    use analysis

    use WRITE_OUT

    implicit none
    real(8) :: start_time, end_time
    real(8) :: time_stamp_1, tmp_time, total_main, total_wall, total_write, total_NNPS
    real(8) :: total_MP, total_ini, total_PS
    !=======================
    ! 1. setting
    !=======================
    !$ start_time = omp_get_wtime()
    write_count = 0
    t = 0
    total_main  = 0.0d0
    total_wall  = 0.0d0
    total_write = 0.0d0
    total_NNPS  = 0.0d0
    total_MP    = 0.0d0
    total_ini   = 0.0d0
    total_PS    = 0.0d0
    call system_setting          ! GLB_SETTING module
    call setting_arrangement     ! INITIAL_SETTING module
    call set_ref_field           ! INITIAL_SETTING module
    call NNPS_ini                ! NNPS module
    call inner_cell_num          ! NNPS module
    call setting_VM              ! INITIAL_SETTING module
    call NNPS_VM                 ! NNPS module
    call NNPS_set_9              ! NNPS module
    call set_MP                  ! INITIAL SETTING module
    call NNPS_MP                 ! NNPS module
    call cal_W_ave               ! cal_KERNEL module
    call write_glb_setting       ! WRITE_OUT module
    call write_initial_setting   ! WRITE_OUT module
    call cal_V_rms_ave           ! ANALYSIS module
    call cal_VM_WALL             ! cal_WALL module
    call write_save_file         ! WRITE_OUT module
    !$ time_stamp_1 = omp_get_wtime()
    total_ini = time_stamp_1 - start_time

    !=======================
    ! 2. time loop
    !=======================
    do t = 1, total_step
        !----------------------------
        ! 2nd-order Adams-Bashforth
        !----------------------------
        !$ time_stamp_1 = omp_get_wtime()
        call cal_main_2AB        ! cal_MAIN module
        !$ tmp_time = omp_get_wtime()
        total_main = total_main + (tmp_time - time_stamp_1)

        ! PS scheme
        if (PS_use == 1) then
            call NNPS_exe        ! NNPS module
            call cal_VM_WALL     ! cal_WALL module

            !$ time_stamp_1 = omp_get_wtime()
            call PS_scheme_exe   ! PS_scheme module
            !$ tmp_time = omp_get_wtime()
            total_PS = total_PS + (tmp_time - time_stamp_1)
        endif

        !$ time_stamp_1 = omp_get_wtime()
        call NNPS_exe            ! NNPS module
         !$ tmp_time = omp_get_wtime()
        total_NNPS = total_NNPS + (tmp_time - time_stamp_1)

        !$ time_stamp_1 = omp_get_wtime()
        call cal_VM_WALL         ! cal_WALL module
        !$ tmp_time = omp_get_wtime()
        total_wall = total_wall + (tmp_time - time_stamp_1)

        !=========================
        ! 3. save data
        !=========================
        if (mod(t, write_step) == 0) then

            !$ time_stamp_1 = omp_get_wtime()
            call interpolate_MP     ! cal_MP module
            !$ tmp_time = omp_get_wtime()
            total_MP = total_MP + (tmp_time - time_stamp_1)

            call cal_V_rms_ave      ! ANALYSIS module

            !$ time_stamp_1 = omp_get_wtime()
            call write_save_file    ! WRITE_OUT module
            !$ tmp_time = omp_get_wtime()
            total_write = total_write + (tmp_time - time_stamp_1)

        endif

        call write_progress(t)
    enddo

    !$ end_time = omp_get_wtime()
    call write_announce_end(start_time, end_time) ! WRITE_OUT module
    call deallocate_exe                           ! GLB_ARG module
    write(*,*) 'time_stamp info :'
    write(*,*) 'initial setting :', total_ini
    write(*,*) 'cal_main        :', total_main
    write(*,*) 'cal_wall        :', total_wall
    write(*,*) 'cal_NNPS        :', total_NNPS
    write(*,*) 'cal_MP          :', total_MP
    write(*,*) 'cal_PS          :', total_PS
    write(*,*) 'writing         :', total_write

end program main