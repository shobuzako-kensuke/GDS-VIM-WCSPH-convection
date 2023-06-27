!==========================================================================================!
!                           thermal convection using GDS-WCSPH                             !
!------------------------------------------------------------------------------------------!
!                          Copyright by Kensuke Shobuzako (2022)                           !
!==========================================================================================!

module GLB_SETTING

    !$ use omp_lib
    
    use INPUT

    implicit none

    ! universal
    real(8), parameter :: pi = dacos(-1.d0)
    integer, save :: i, j, k, t, tmp_int, me, you
    real(8), save :: tmp

    ! saved
    real(8), save :: total_num, length, T_bot, T_top, rho_0, alpha, c_p, k_h, eta_0, K_0
    real(8), save :: delta_x, smooth_len, SS, RSS, delta_tem, T_0
    real(8), save :: SP_mass, kappa, nu, Ra, Pr, M_th
    real(8), save :: t_CFL, r_t_CFL, t_VNM, r_t_VNM, t_VNT, time_step
    real(8), save :: alpha_d, sprt_dom
    real(8), save :: pass_time
    integer, save :: accel_step_end

    ! tmp
    real(8) :: system_V, eff_Pr
    real(8) :: V_fin, V_fin_tb
    real(8) :: tau_V_fin, tau_th, step_V_fin, step_th, eff_write, total_time
    
contains

    !---------------------
    ! system setting
    !---------------------
    subroutine system_setting
        if (RSST_VIM == 0) then
            length = length_R
            T_top  = T_top_R
            if (Ra_0 == 1.0d+4) then
                delta_tem = 1.0d-2
            elseif (Ra_0 == 1.0d+5) then
                delta_tem = 1.0d-1
            elseif(Ra_0 == 1.0d+6) then
                delta_tem = 1.0d0
            endif
            T_bot  = T_top + delta_tem
            rho_0  = rho_0_R
            k_h    = k_h_R
            eta_0  = eta_R
            alpha  = alpha_R
            c_p    = c_p_R
            K_0    = K_0_R
        elseif (RSST_VIM == 1) then
            length = length_V
            T_top  = T_top_V
            T_bot  = T_bot_V
            delta_tem = T_bot - T_top
            if (Ra_0 == 1.0d+4) then
                eta_0 = 1.0d+23
            elseif (Ra_0 == 1.0d+5) then
                eta_0 = 1.0d+22
            elseif (Ra_0 == 1.0d+6) then
                eta_0 = 1.0d+21
            endif
            rho_0  = rho_0_V
            k_h    = k_h_V
            alpha  = alpha_V
            c_p    = c_p_V
            K_0    = K_0_V
        endif

        total_num  = num**2                       ! the number of total SP
        delta_x    = length / dble(num)           ! regular particle distance
        smooth_len = coe_smooth * delta_x         ! smoothing length
        SS         = (K_0/rho_0) ** (1.0d0/2.0d0) ! speed of sound [m s-2]
        RSS        = SS / (zeta * xi)             ! reduced SS [m s-2]
        T_0        = 0.5d0 * (T_bot + T_top)      ! reference temperature [K]
        system_V  = length**2                     ! system volume
        
        SP_mass   = (rho_0 * system_V) / dble(total_num)       ! particle mass
        nu        = eta_0 / rho_0                              ! kinematic viscosity [m2 s-1]
        kappa     = k_h / (rho_0*c_p)                          ! thermal diffusivity [m2 s-1]
        Pr        = nu / kappa                                 ! Prandtl number
        eff_Pr    = Pr / (xi**2)                               ! effective Prandtl number
        Ra        = (alpha * (-g) * delta_tem * (length)**3) & ! Rayleigh number
                    / (kappa * nu)
        M_th      = kappa / (length * SS)                      ! Mach number
        t_CFL     = coe_step*smooth_len / SS                   ! CFL condition [s]
        r_t_CFL   = coe_step*smooth_len / RSS                  ! reduced CFL condition [s]
        t_VNM     = coe_step*(smooth_len**2) / (2.0d0 * nu)    ! von Neumann condition momentum
        r_t_VNM   = t_VNM * xi**2                              ! reduced von Neumann condition for EOM
        t_VNT     = coe_step*(smooth_len**2) / (2.0d0 * kappa) ! von Neumann condition for thermal diffusion

        ! time step
        time_step = min(r_t_CFL, r_t_VNM, t_VNT)
        
        if (RSST_VIM == 0) then
            ! free-fall velocity
            v_fin    = (Ra * Pr)**(1.0d0/2.0d0) * M_th * SS
            v_fin_tb = Ra**(3.0d0/7.0d0) * Pr**(2.0d0/7.0d0) * M_th * SS  ! new scaling(2023/06/22)
            !v_fin_tb = (Ra * Pr)**(2.0d0/5.0d0) * M_th * SS ! for syuron scaling

        elseif (RSST_VIM == 1) then
            ! stokes velocity
            v_fin    = Ra * M_th * SS
            v_fin_tb = Ra**(2.0d0/3.0d0) * M_th * SS
        endif

        tau_V_fin  = length / V_fin                ! typical time scale for terminal velocity [s]       
        tau_th     = length**2 / kappa             ! typical time scale for thermal diffusion
        step_V_fin = tau_V_fin / time_step         ! typical steps for terminal velocity
        step_th    = tau_th / time_step            ! typical steps for thermal diffusion
        eff_write  = dble(write_step) * time_step  ! effective writing interval [s] 
        total_time = time_step * dble(total_step)  ! total time [s]
        
        if (kernel_switch == 0) then
            alpha_d  = 7.0d0 / (478.0d0 * pi * smooth_len**2) ! Quintic
            sprt_dom = 3.0d0 * smooth_len
        elseif (kernel_switch == 1) then
            alpha_d  = 9.0d0 / (4.0d0   * pi * smooth_len**2) ! Wendland C4
            sprt_dom = 2.0d0 * smooth_len
        elseif (kernel_switch == 2) then
            alpha_d  = 39.0d0 / (14.0d0 * pi * smooth_len**2) ! Wendland C6
            sprt_dom = 2.0d0 * smooth_len
        endif
    end subroutine system_setting
    
end module GLB_SETTING