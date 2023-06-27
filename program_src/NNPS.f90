!==========================================================================================!
!                           thermal convection using GDS-WCSPH                             !
!------------------------------------------------------------------------------------------!
!                          Copyright by Kensuke Shobuzako (2022)                           !
!==========================================================================================!

module NNPS
    
    use INPUT
    use GLB_SETTING
    use GLB_ARG
    use INITIAL_SETTING

    implicit none

    integer, save :: cell_cap
    integer, save :: inner_cell_start, inner_cell_end
    integer :: cell_num

contains
    !---------------------
    ! NNPS setting
    !---------------------
    subroutine NNPS_ini
        ! max cell number
        NNPS_max = ceiling(SP_x_max/sprt_dom) + (ceiling(SP_y_max/sprt_dom) - 1) &
                                               * ceiling(SP_x_max/sprt_dom)
        ! set NNPS_arg
        tmp = sprt_dom / delta_x
        cell_cap = (ceiling(tmp) + 2)**2
        allocate(NNPS_arg(cell_cap, NNPS_max))
        NNPS_arg(:,:) = 0

        ! set NNPS_row
        allocate(NNPS_row(NNPS_max))
        NNPS_row(:) = 1

        do i = 1, total_num_system
            ! belonging cell number
            cell_num = ceiling(SP_xy(i, 1)/sprt_dom) + (ceiling(SP_xy(i, 2)/sprt_dom) - 1) &
                                                     * ceiling(SP_x_max   /sprt_dom)
            ! into NNPS_arg
            NNPS_arg(NNPS_row(cell_num), cell_num) = i ! substitute particle number
            NNPS_row(cell_num) = NNPS_row(cell_num) + 1

            ! check program
            if (NNPS_row(cell_num) > cell_cap) stop &
                '[Error] over NNPS capacity at initial NNPS'

            ! save for initial state
            SP_uv(i, 1) = cell_num
        enddo
    end subroutine NNPS_ini

    !---------------------
    ! inner cell num
    !---------------------
    subroutine inner_cell_num
        inner_cell_start = 0
        inner_cell_end   = 0
        do i = 1, total_num_system
            if ((SP_kind(i) == 0) .and. (inner_cell_start == 0)) then
                inner_cell_start = ceiling(SP_xy(i, 1)/sprt_dom) + (ceiling(SP_xy(i, 2)/sprt_dom) - 1) &
                                                                  * ceiling(SP_x_max   /sprt_dom)
            elseif (SP_kind(i) == 0) then
                tmp_int = ceiling(SP_xy(i, 1)/sprt_dom) + (ceiling(SP_xy(i, 2)/sprt_dom) - 1) &
                                                         * ceiling(SP_x_max   /sprt_dom)
                if (tmp_int < inner_cell_start) then
                    inner_cell_start = tmp_int
                elseif (tmp_int > inner_cell_end) then
                    inner_cell_end = tmp_int
                endif
            endif
        enddo
    end subroutine inner_cell_num

    !---------------------
    ! NNPS for VM
    !---------------------
    subroutine NNPS_VM
        ! set NNPS_arg_VM
        allocate(NNPS_arg_VM(cell_cap*2, NNPS_max))
        NNPS_arg_VM(:,:) = 0

        ! set NNPS_VM_row
        allocate(NNPS_VM_row(NNPS_max))
        NNPS_VM_row(:) = 1

        do i = 1, num_VM
            ! belonging cell number
            cell_num = ceiling(VM_xy(i, 1)/sprt_dom) + (ceiling(VM_xy(i, 2)/sprt_dom) - 1) &
                                                     * ceiling(SP_x_max   /sprt_dom)
            ! into NNPS_arg_VM
            NNPS_arg_VM(NNPS_VM_row(cell_num), cell_num) = i ! substitute particle number                 
            NNPS_VM_row(cell_num) = NNPS_VM_row(cell_num) + 1

            ! checl program
            if (NNPS_VM_row(cell_num) > (cell_cap*2)) stop &
                '[Error] over NNPS capacity at initial NNPS for VM'
            
            ! save for initial state
            VM_uv(i, 1) = cell_num
        enddo
    end subroutine NNPS_VM

    !---------------------
    ! set 9 cells
    !---------------------
    subroutine NNPS_set_9
        allocate(NNPS_9(9, NNPS_max))
        NNPS_9(:,:) = 0
        !$ call omp_set_dynamic(.false.)
        !$ call omp_set_num_threads(num_threads)
        !$OMP parallel default(none) &
        !$OMP private(i) &
        !$OMP shared(NNPS_max, NNPS_9, SP_x_max, sprt_dom)
        !$OMP do
        do i = 1, NNPS_max
            ! NNPS_9(9 cells, my cell number)
            NNPS_9(5, i) = i
            NNPS_9(4, i) = NNPS_9(5, i) - 1
            NNPS_9(6, i) = NNPS_9(5, i) + 1
            NNPS_9(2, i) = NNPS_9(5, i) - ceiling(SP_x_max/sprt_dom)
            NNPS_9(1, i) = NNPS_9(2, i) - 1
            NNPS_9(3, i) = NNPS_9(2, i) + 1
            NNPS_9(8, i) = NNPS_9(5, i) + ceiling(SP_x_max/sprt_dom)
            NNPS_9(7, i) = NNPS_9(8, i) - 1
            NNPS_9(9, i) = NNPS_9(8, i) + 1 
        enddo
        !$OMP enddo
        !$OMP end parallel
    end subroutine NNPS_set_9

    !---------------------
    ! NNPS_mapping
    !---------------------
    subroutine NNPS_MP
        ! set NNPS_arg_MP
        allocate(NNPS_arg_MP(cell_cap, NNPS_max))
        NNPS_arg_MP(:, :) = 0

        ! into NNPS_arg_MP
        do i = 1, num_MP
            cell_num = ceiling(MP_xy(i, 1)/sprt_dom) + (ceiling(MP_xy(i, 2)/sprt_dom) - 1) &
                                                     * ceiling(SP_x_max   /sprt_dom)
            do j = 1, cell_cap
                if (NNPS_arg_MP(j, cell_num) == 0) then
                    NNPS_arg_MP(j, cell_num) = i ! substitute particle number
                    exit
                endif

                ! check program
                if (j > cell_cap) stop &
                        '[Error] over NNPS capacity at initial NNPS for MP'
            enddo
            
            ! save for initial state
            MP_uv(i, 1) = cell_num
        enddo
    end subroutine NNPS_MP

    !---------------------
    ! NNPS in time_loop
    !---------------------
    subroutine NNPS_exe
        ! clear
        NNPS_arg(:,:) = 0
        NNPS_row(:)   = 1
        do i = 1, total_num_system
            ! belonging cell number
            cell_num = ceiling(SP_xy(i, 1)/sprt_dom) + (ceiling(SP_xy(i, 2)/sprt_dom) - 1) &
                                                     * ceiling(SP_x_max   /sprt_dom)
            
            ! ! check program
            ! if (NNPS_row(cell_num) > (cell_cap)) stop '[Error] over NNPS capacity'

            ! into NNPS_arg
            NNPS_arg(NNPS_row(cell_num), cell_num) = i ! substitute particle number
            NNPS_row(cell_num) = NNPS_row(cell_num) + 1
        enddo
    end subroutine NNPS_exe

end module NNPS