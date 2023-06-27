!==========================================================================================!
!                           thermal convection using GDS-WCSPH                             !
!------------------------------------------------------------------------------------------!
!                          Copyright by Kensuke Shobuzako (2022)                           !
!==========================================================================================!

module GLB_ARG
    use INPUT
    implicit none

    ! particles
    real(8), allocatable, save :: SP_xy(:,:)       ! position
    real(8), allocatable, save :: SP_uv(:,:)       ! velocity
    real(8), allocatable, save :: SP_rho_S(:)      ! SP density
    real(8), allocatable, save :: SP_pre(:)        ! pressure
    real(8), allocatable, save :: SP_tem(:)        ! temperature
    real(8), allocatable, save :: SP_eta(:)        ! viscosity
    real(8), allocatable, save :: SP_sm_rho(:)     ! smoothed density
    real(8), allocatable, save :: SP_rho_G(:)      ! gravitational density
    real(8), allocatable, save :: par_con(:)       ! particle consistency

    ! classify particles
    integer, allocatable, save :: SP_kind(:)       ! 0(inner), 1-4(wall)

    ! NNPS
    integer, allocatable, save :: NNPS_arg(:,:)    ! (SP_num, cell_num)
    integer, allocatable, save :: NNPS_9(:,:)      ! (9 neighber_cell_nums, cell_num)
    integer, allocatable, save :: NNPS_row(:)      ! (the number of particles in each cell)

    ! virtual markers
    real(8), allocatable, save :: VM_xy(:,:)       ! position
    real(8), allocatable, save :: VM_uv(:,:)       ! velocity
    real(8), allocatable, save :: VM_pre(:)        ! pressure
    real(8), allocatable, save :: VM_tem(:)        ! temperature
    real(8), allocatable, save :: VM_rho(:)        ! density
    real(8), allocatable, save :: VM_corrected(:)  ! denominator for correcting VM
    real(8), allocatable, save :: VM_F_i(:,:)      ! F_i 

    ! ID between virtual markers and wall particles
    integer, allocatable, save :: ID_VM_WL(:)      ! your wall number
    integer, allocatable, save :: NNPS_arg_VM(:,:) ! (VM_num, cell_num)
    integer, allocatable, save :: NNPS_VM_row(:)   ! (the number of VM in each cell)

    ! EOC, EOM, EOE
    real(8), allocatable, save :: EOC_b(:)         ! before EOC
    real(8), allocatable, save :: EOC_n(:)         ! now EOC
    real(8), allocatable, save :: EOM_b(:,:)       ! before EOM
    real(8), allocatable, save :: EOM_n(:,:)       ! now EOM
    real(8), allocatable, save :: EOE_b(:)         ! before EOE
    real(8), allocatable, save :: EOE_n(:)         ! now EOE
    real(8), allocatable, save :: vis_EOM(:,:)     ! viscosity term in EOM
    real(8), allocatable, save :: buo_EOM(:)       ! buoyancy term in EOM
    real(8), allocatable, save :: pre_EOM(:,:)     ! pressure gradient term in EOM
    
    ! F_i for VM
    real(8), allocatable, save :: F_i(:,:)         ! inner particle EOM for VM

    ! PS scheme
    real(8), allocatable, save :: PS_shift_xy(:,:) ! shift vector
    real(8), allocatable, save :: PS_TE_uv(:,:)    ! using Taylor expansion 
    real(8), allocatable, save :: PS_TE_rho_S(:,:)  
    real(8), allocatable, save :: PS_TE_pre(:,:)
    real(8), allocatable, save :: PS_TE_tem(:,:)
    real(8), allocatable, save :: PS_TE_sm_rho(:,:)
    real(8), allocatable, save :: PS_KGC(:,:)      ! KGC arrangement

    ! mapping
    real(8), allocatable, save :: MP_xy(:,:)       ! position
    real(8), allocatable, save :: MP_uv(:,:)       ! velocity
    real(8), allocatable, save :: MP_pre(:)        ! pressure
    real(8), allocatable, save :: MP_tem(:)        ! temperature
    real(8), allocatable, save :: MP_NU(:)         ! Nu number
    integer, allocatable, save :: NNPS_arg_MP(:,:) ! NNPS for MP

    ! mapping KGC
    real(8), allocatable, save :: MP_KGC(:,:)      ! KGC arrangement
    real(8), allocatable, save :: MP_corrected(:)  ! denominator for correcting MP

contains
    subroutine deallocate_exe
        deallocate(SP_xy, SP_uv, SP_rho_S, SP_pre, SP_tem, SP_eta, SP_sm_rho, SP_rho_G, par_con)
        deallocate(SP_kind, NNPS_arg, NNPS_9, NNPS_row)
        deallocate(VM_xy, VM_uv, VM_pre, VM_tem, VM_rho, VM_corrected, VM_F_i)
        deallocate(ID_VM_WL, NNPS_arg_VM, NNPS_VM_row)
        deallocate(EOC_b, EOC_n, EOM_b, EOM_n, EOE_b, EOE_n, vis_EOM, buo_EOM, pre_EOM)
        deallocate(F_i)
        deallocate(MP_xy, MP_uv, MP_pre, MP_tem, NNPS_arg_MP)
        deallocate(MP_KGC, MP_corrected, MP_NU)
        deallocate(PS_shift_xy, PS_TE_uv, PS_TE_rho_S, PS_TE_pre, PS_TE_tem, PS_TE_sm_rho, PS_KGC)
        
    end subroutine deallocate_exe

end module GLB_ARG