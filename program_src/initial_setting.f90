!==========================================================================================!
!                           thermal convection using GDS-WCSPH                             !
!------------------------------------------------------------------------------------------!
!                          Copyright by Kensuke Shobuzako (2022)                           !
!==========================================================================================!

module INITIAL_SETTING

    !$ use omp_lib
    use INPUT
    use GLB_SETTING
    use GLB_ARG
    use cal_EOS

    implicit none
    integer, save :: wall_thickness, num_system, total_num_system
    real(8), save :: length_system, SP_x_max, SP_y_max, wall_width
    integer, save :: num_VM
    integer, save :: NNPS_max
    integer, parameter :: num_MP = num_MP_side*6
    real(8) :: MP_delta_x

contains

    !---------------------
    ! set position
    !---------------------
    subroutine setting_arrangement
        ! define wall thickness
        wall_thickness = ceiling(sprt_dom/delta_x) + 2
        ! the num of inner and wall particles in an edge
        num_system = num + 2 * wall_thickness
        ! the total_num of inner and wall particles in the system
        total_num_system = num_system**2
        ! system length
        length_system = length + 2.0d0 * (dble(wall_thickness) * delta_x)

        ! allocate
        allocate(SP_xy    (total_num_system, 2))
        allocate(SP_uv    (total_num_system, 2))
        allocate(SP_rho_S (total_num_system))
        allocate(SP_pre   (total_num_system))
        allocate(SP_tem   (total_num_system))
        allocate(SP_eta   (total_num_system))
        allocate(SP_sm_rho(total_num_system))
        allocate(SP_rho_G (total_num_system))
        allocate(SP_kind  (total_num_system))
        allocate(par_con  (total_num_system))
        allocate(vis_EOM  (total_num_system, 2))
        allocate(buo_EOM  (total_num_system))
        allocate(pre_EOM  (total_num_system, 2))
    
        allocate(EOC_b(total_num_system))
        allocate(EOC_n(total_num_system))
        allocate(EOM_b(total_num_system, 2))
        allocate(EOM_n(total_num_system, 2))
        allocate(EOE_b(total_num_system))
        allocate(EOE_n(total_num_system))
        allocate(F_i  (total_num_system, 2))

        ! for Taylor expansion in PS scheme
        allocate(PS_shift_xy (total_num_system, 2))
        allocate(PS_TE_uv    (total_num_system, 4))
        allocate(PS_TE_rho_S (total_num_system, 2))
        allocate(PS_TE_pre   (total_num_system, 2))
        allocate(PS_TE_tem   (total_num_system, 2))
        allocate(PS_TE_sm_rho(total_num_system, 2))
        allocate(PS_KGC      (total_num_system, 4))

        SP_xy(:,:)   = 0.0d0
        SP_uv(:,:)   = 0.0d0
        SP_rho_S(:)  = rho_0
        SP_pre(:)    = 0.0d0
        SP_tem(:)    = 0.0d0
        SP_eta(:)    = eta_0
        SP_sm_rho(:) = 0.0d0
        SP_rho_G(:)  = rho_0
        SP_kind(:)   = 1      ! if all particles were wall
        par_con(:)   = 0.0d0
        vis_EOM(:,:) = 0.0d0
        buo_EOM(:)   = 0.0d0
        pre_EOM(:,:) = 0.0d0
        
        PS_shift_xy(:,:) = 0.0d0
        
        EOC_b(:)   = 0.0d0
        EOC_n(:)   = 0.0d0
        EOM_b(:,:) = 0.0d0
        EOM_n(:,:) = 0.0d0
        EOE_b(:)   = 0.0d0
        EOE_n(:)   = 0.0d0

        ! initial regular arrangement
        SP_xy(1, 1) = delta_x / 2.0d0
        SP_xy(1, 2) = delta_x / 2.0d0
        do i = 2, total_num_system
            if (i <= num_system) then
                SP_xy(i, 1) = SP_xy(i-1, 1) + delta_x
                SP_xy(i, 2) = SP_xy(i-1, 2)
            else
                SP_xy(i, 1) = SP_xy(i-num_system, 1)
                SP_xy(i, 2) = SP_xy(i-num_system, 2) + delta_x
            endif
        enddo
        SP_x_max = maxval(SP_xy(:, 1))
        SP_y_max = maxval(SP_xy(:, 2))

        ! classify particles
        wall_width = dble(wall_thickness) * delta_x
        !$ call omp_set_dynamic(.false.)
        !$ call omp_set_num_threads(num_threads)
        !$OMP parallel default(shared) &
        !$OMP private(i) & 
        !$OMP shared(total_num_system, SP_xy, length, SP_kind, wall_width)
        !$OMP do
        do i = 1, total_num_system
            ! find inner particles
            if ((wall_width < SP_xy(i, 1)) .and. (SP_xy(i, 1) < length + wall_width) .and. &
                (wall_width < SP_xy(i, 2)) .and. (SP_xy(i, 2) < length + wall_width)) then
                SP_kind(i) = 0
            endif
            
            ! classify wall particles
            if (SP_kind(i) /= 0) then
                if (SP_xy(i, 2) < wall_width) then
                    SP_kind(i) = 1 ! bottom
                elseif (SP_xy(i, 2) > length + wall_width) then
                    SP_kind(i) = 2 ! top
                elseif (SP_xy(i, 1) > length + wall_width) then
                    SP_kind(i) = 3 ! right
                else
                    SP_kind(i) = 4 ! left                        
                endif
            endif
        enddo
        !$OMP end do
        !$OMP barrier
        !$OMP end parallel
    end subroutine setting_arrangement

    !---------------------
    ! set reference field
    !---------------------
    subroutine set_ref_field
        if ((RSST_VIM == 0) .or. (RSST_VIM == 1)) then
            ! linear temperature profile
            SP_tem(:) = T_bot - (delta_tem / length) * (SP_xy(:, 2) - wall_width)

            ! temperature perturbation
            !$ call omp_set_dynamic(.false.)
            !$ call omp_set_num_threads(num_threads)
            !$OMP parallel default(none) &
            !$OMP private(i) &
            !$OMP shared(total_num_system, SP_xy, length, wall_width, SP_tem, delta_tem) &
            !$OMP shared(SP_kind, T_bot, T_top, SP_pre, SP_rho_G)
            !$OMP do
            do i = 1, total_num_system
                ! set temperature perturbation
                if ((SP_xy(i, 2) - wall_width) / length < 0.1) then
                    SP_tem(i) = SP_tem(i) + 0.3d0 * delta_tem * &
                                dcos(pi * (SP_xy(i, 1) - wall_width) / length)
                endif

                ! if bot or top wall, fixed temperature
                if (SP_kind(i) == 1) then
                    SP_tem(i) = T_bot
                elseif (SP_kind(i) == 2) then
                    SP_tem(i) = T_top
                endif

                ! cal rho_G
                call cal_real_EOS(SP_pre(i), SP_tem(i), SP_rho_G(i))
            enddo
            !$OMP enddo
            !$OMP barrier
            !$OMP end parallel
        endif
    end subroutine set_ref_field    

    !---------------------
    ! set virtual markers
    !---------------------
    subroutine setting_VM
        ! define the number of VM
        tmp_int = 0
        do i = 1, total_num_system
            if (SP_kind(i) /= 0) then
                tmp_int = tmp_int + 1
            endif
        enddo
        
        ! check program
        if (tmp_int /= (total_num_system - total_num)) stop &
            '[Error] num_VM not equal num_wall at initial_setting'

        ! allocate
        num_VM = tmp_int
        allocate(VM_xy       (num_VM, 2))
        allocate(VM_uv       (num_VM, 2))
        allocate(VM_pre      (num_VM))
        allocate(VM_tem      (num_VM))
        allocate(VM_rho      (num_VM))
        allocate(ID_VM_WL    (num_VM))
        allocate(VM_corrected(num_VM))
        allocate(VM_F_i      (num_VM, 2))
        VM_xy(:,:)  = 0.0d0
        VM_uv(:,:)  = 0.0d0
        VM_pre(:)   = 0.0d0
        VM_tem(:)   = 0.0d0
        VM_rho(:)   = 0.0d0
        ID_VM_WL(:) = 0

        ! make pairs
        j = 1
        do i = 1, total_num_system
            if (SP_kind(i) /= 0) then
                ID_VM_WL(j) = i
                j = j + 1

                ! check program
                if (j > (num_VM+1)) stop &
                    '[Error] over ID_VM_WL at initial_setting'
            endif
        enddo

        ! position
        !$ call omp_set_dynamic(.false.)
        !$ call omp_set_num_threads(num_threads)
        !$OMP parallel default(none) &
        !$OMP private(me, you) &
        !$OMP shared(num_VM, ID_VM_WL, SP_kind, SP_xy, VM_xy, length, wall_width)
        !$OMP do
        do me = 1, num_VM
            you = ID_VM_WL(me) ! your wall number
            ! bottom wall
            if (SP_kind(you) == 1) then
                if (SP_xy(you, 1) < wall_width) then
                    VM_xy(me, 1) = 2.0d0 * wall_width - SP_xy(you, 1)
                    VM_xy(me, 2) = 2.0d0 * wall_width - SP_xy(you, 2)
                elseif ((length + wall_width) < SP_xy(you, 1)) then
                    VM_xy(me, 1) = 2.0d0 * (wall_width + length) - SP_xy(you, 1)
                    VM_xy(me, 2) = 2.0d0 * (wall_width         ) - SP_xy(you, 2)
                else
                    VM_xy(me, 1) = SP_xy(you, 1)
                    VM_xy(me, 2) = 2.0d0 * wall_width - SP_xy(you, 2)
                endif
            ! top wall
            elseif (SP_kind(you) == 2) then
                if (SP_xy(you, 1) < wall_width) then
                    VM_xy(me, 1) = 2.0d0 * (wall_width         ) - SP_xy(you, 1)
                    VM_xy(me, 2) = 2.0d0 * (wall_width + length) - SP_xy(you, 2)
                elseif ((length + wall_width) < SP_xy(you, 1)) then
                    VM_xy(me, 1) = 2.0d0 * (wall_width + length) - SP_xy(you, 1)
                    VM_xy(me, 2) = 2.0d0 * (wall_width + length) - SP_xy(you, 2)
                else
                    VM_xy(me, 1) = SP_xy(you, 1)
                    VM_xy(me, 2) = 2.0d0 * (length + wall_width) - SP_xy(you, 2)
                endif
            ! right wall
            elseif(SP_kind(you) == 3) then
                VM_xy(me, 1) = 2.0d0 * (length + wall_width) - SP_xy(you, 1)
                VM_xy(me, 2) = SP_xy(you, 2)
            else
                VM_xy(me, 1) = 2.0d0 * wall_width - SP_xy(you, 1)
                VM_xy(me, 2) = SP_xy(you, 2)

            ! check program
            if (SP_kind(you) == 0) stop &
                '[Error] Inner particle was selected to VM'

            endif
        enddo
        !$OMP enddo
        !$OMP end parallel
    end subroutine setting_VM

    !---------------------
    ! set mapping
    !---------------------
    subroutine set_MP
        ! set 50 calculation points in each edge
        allocate(MP_xy       (num_MP, 2))
        allocate(MP_uv       (num_MP, 2))
        allocate(MP_pre      (num_MP))
        allocate(MP_tem      (num_MP))
        allocate(MP_KGC      (num_MP, 4))
        allocate(MP_corrected(num_MP))
        allocate(MP_NU       (num_MP))
        MP_xy(:,:) = 0.0d0
        MP_uv(:,:) = 0.0d0
        MP_pre(:)  = 0.0d0
        MP_tem(:)  = 0.0d0
        MP_NU(:)   = 0.0d0

        !---------------------
        ! position & kind
        !---------------------
        MP_delta_x = length / (dble(num_MP_side-1))

        do i = 1, num_MP
            ! bottom layer
            if ((0*num_MP_side+1 <= i) .and. (i < 1*num_MP_side+1)) then
                if (i == 0*num_MP_side+1) then
                    MP_xy(i, 1) = wall_width
                    MP_xy(i, 2) = wall_width
                else
                    MP_xy(i, 1) = MP_xy(i-1, 1) + MP_delta_x
                    MP_xy(i, 2) = MP_xy(i-1, 2)
                endif
            ! horizontal middle layer
            elseif ((1*num_MP_side+1 <= i) .and. (i < 2*num_MP_side+1)) then
                if (i == 1*num_MP_side+1) then
                    MP_xy(i, 1) = wall_width
                    MP_xy(i, 2) = wall_width + length / 2.0d0
                else
                    MP_xy(i, 1) = MP_xy(i-1, 1) + MP_delta_x
                    MP_xy(i, 2) = MP_xy(i-1, 2)
                endif
            ! top layer
            elseif ((2*num_MP_side+1 <= i) .and. (i < 3*num_MP_side+1)) then
                if (i == 2*num_MP_side+1) then
                    MP_xy(i, 1) = wall_width
                    MP_xy(i, 2) = wall_width + length
                else
                    MP_xy(i, 1) = MP_xy(i-1, 1) + MP_delta_x
                    MP_xy(i, 2) = MP_xy(i-1, 2)
                endif
            ! left layer
            elseif ((3*num_MP_side+1 <= i) .and. (i < 4*num_MP_side+1)) then
                if (i == 3*num_MP_side+1) then
                    MP_xy(i, 1) = wall_width
                    MP_xy(i, 2) = wall_width
                else
                    MP_xy(i, 1) = MP_xy(i-1, 1)
                    MP_xy(i, 2) = MP_xy(i-1, 2) + MP_delta_x
                endif
            ! vertical middle layer
            elseif ((4*num_MP_side+1 <= i) .and. (i < 5*num_MP_side+1)) then
                if (i == 4*num_MP_side+1) then
                    MP_xy(i, 1) = wall_width + length / 2.0d0
                    MP_xy(i, 2) = wall_width 
                else
                    MP_xy(i, 1) = MP_xy(i-1, 1)
                    MP_xy(i, 2) = MP_xy(i-1, 2) + MP_delta_x
                endif
            ! right layer
            elseif ((5*num_MP_side+1 <= i) .and. (i < 6*num_MP_side+1)) then
                if (i == 5*num_MP_side+1) then
                    MP_xy(i, 1) = wall_width + length
                    MP_xy(i, 2) = wall_width 
                else
                    MP_xy(i, 1) = MP_xy(i-1, 1)
                    MP_xy(i, 2) = MP_xy(i-1, 2) + MP_delta_x
                endif
            endif
        enddo
    end subroutine set_MP

end module INITIAL_SETTING