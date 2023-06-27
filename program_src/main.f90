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
    use set_time_step
    use cal_WALL
    use cal_MAIN
    use PS_scheme
    use cal_MAPPING
    use analysis

    use WRITE_OUT

    implicit none
    real(8) :: start_time, end_time
    !=======================
    ! 1. setting
    !=======================
    !$ start_time = omp_get_wtime()
    write_count = 0
    t = 0
    pass_time = 0.0d0
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
    call analysis_exe            ! ANALYSIS module
    call cal_VM_WALL             ! cal_WALL module
    call write_save_file         ! WRITE_OUT module
    
    !=======================
    ! 2. time loop
    !=======================
    do t = 1, total_step

        call set_delta_t         ! set_time_step module
        !----------------------------
        ! calculation EOM and EOE
        ! 2nd-order Adams-Bashforth
        !----------------------------
        call cal_main_2AB        ! cal_MAIN module
        !----------------------------
        ! PS Technology
        !----------------------------
        call NNPS_exe            ! NNPS module
        call cal_VM_WALL         ! cal_WALL module
        call PS_scheme_exe       ! PS_scheme module
        !----------------------------
        ! renewal wall
        !----------------------------
        call NNPS_exe            ! NNPS module
        call cal_VM_WALL         ! cal_WALL module

        !=========================
        ! 3. save data
        !=========================
        if (mod(t, write_step) == 0) then
            call interpolate_MP     ! cal_MP module
            call analysis_exe       ! ANALYSIS module
            call write_save_file    ! WRITE_OUT module
        endif

        call write_progress(t)      ! WRITE_OUT module
    enddo

    !$ end_time = omp_get_wtime()
    call write_announce_end(start_time, end_time) ! WRITE_OUT module
    call deallocate_exe                           ! GLB_ARG module

end program main