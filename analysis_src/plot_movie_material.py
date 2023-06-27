#==========================================================================================#
#                           thermal convection using GDS-WCSPH                             #
#------------------------------------------------------------------------------------------#
#                          Copyright by Kensuke Shobuzako (2022)                           #
#==========================================================================================#

#---------------------
# CHARM
#---------------------
import numpy as np
import matplotlib.pyplot as plt
import math
import glob
from scipy import integrate
from scipy.interpolate import interp1d
from matplotlib import rc
import sys
import os
import time
rc('text', usetex=True)

#=====================
# import parameter
#=====================
import read_parameter as rp
import input

plot_without_wall     = input.plot_without_wall
plot_with_wall        = input.plot_with_wall
plot_all_without_wall = input.plot_all_without_wall
plot_all_with_wall    = input.plot_all_with_wall
plot_PS_shift         = input.plot_PS_shift

plot_phy = input.plot_phy

num_save_file    = input.num_save_file
total_step       = rp.total_step
num              = rp.num
num_system       = rp.num_system
total_num        = rp.total_num
total_num_system = rp.total_num_system
num_wall         = rp.num_wall
num_MP_side      = rp.num_MP_side
num_MP           = rp.num_MP
wall_thickness   = rp.wall_thickness
write_step       = rp.write_step
RSST_VIM         = rp.RSST_VIM
kernel_switch    = rp.kernel_switch
 
length         = rp.length
delta_tem      = rp.delta_tem
coe_smooth     = rp.coe_smooth
delta_x        = rp.delta_x
smooth_len     = rp.smooth_len
support_domain = rp.support_domain
time_step      = rp.time_step
 
rho_0   = rp.rho_0
alpha   = rp.alpha
eta_0   = rp.eta_0
T_0     = rp.T_0
T_top   = rp.T_top
kappa   = rp.kappa
 
V_fin    = rp.V_fin
V_fin_tb = rp.V_fin_tb

Ra   = rp.Ra
Pr   = rp.Pr
M_th = rp.M_th
SS   = rp.SS
RSS  = rp.RSS

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
#                          program below                              #
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
print('-----------------------------------------------------')
print('[Message] plot_movie_material has started.')
print('-----------------------------------------------------')
start_time = time.perf_counter()

#---------------------
# scatter size
#---------------------
if (num == 50):
    SP_size = 18
elif (num == 60):
    SP_size = 12
elif (num == 100):
    SP_size = 4
elif (num == 120):
    SP_size = 2
elif (num == 200):
    SP_size = 0.2
elif (num == 240):
    SP_size = 0.13
else:
    SP_size = 10

#---------------------
# comparison set
#---------------------
# normalizing velocity
if (RSST_VIM == 0):
        V_nor = V_fin
elif (RSST_VIM == 1):
        V_nor = kappa / length

# lim settings
if (RSST_VIM == 0):
        if (Ra == 10**4):
                V_max_plot = 0.28
        elif (Ra == 10**5):
                V_max_plot = 0.38
        elif (Ra == 10**6):
                V_max_plot = 0.40
elif (RSST_VIM == 1):
        if (Ra == 10**4):
                V_max_plot = 50
        elif (Ra == 10**5):
                V_max_plot = 200
        elif (Ra == 10**6):
                V_max_plot = 850
# normalizing pressure
pre_nor = rho_0 * alpha * delta_tem * length * 10 
#---------------------
# read data
#---------------------
# get file_name
file_path = glob.glob('./../output_setting/initial_state_SP_*') # finding file
file_name = file_path[0][-8:-4] # [-8:-4] is file_name in str type
endian = '>'

# SP
SP_read = np.zeros((total_num_system, 11, num_save_file+1)) # num_save_file = initial_state + save_state
for i in range(num_save_file+1):
        file_num = str(i)
        f = open('./../save_file/SP_save/SP_save_{}_{}.dat'.format(file_name, file_num), 'rb')
        data_type = np.dtype([('SP_xy',     endian+str(2*total_num_system)+'d'), \
                              ('SP_uv',     endian+str(2*total_num_system)+'d'), \
                              ('SP_rho_S',  endian+str(total_num_system)+'d'), \
                              ('SP_pre',    endian+str(total_num_system)+'d'), \
                              ('SP_tem',    endian+str(total_num_system)+'d'), \
                              ('SP_eta',    endian+str(total_num_system)+'d'), \
                              ('SP_sm_rho', endian+str(total_num_system)+'d'), \
                              ('SP_rho_G',  endian+str(total_num_system)+'d'), \
                              ('par_con',   endian+str(total_num_system)+'d')])
        data = np.fromfile(f, dtype=data_type)
        SP_xy_tmp     = data['SP_xy'].reshape(    (total_num_system, 2), order='F')
        SP_uv_tmp     = data['SP_uv'].reshape(    (total_num_system, 2), order='F')
        SP_rho_S_tmp  = data['SP_rho_S'].reshape( (total_num_system, 1), order='F')
        SP_pre_tmp    = data['SP_pre'].reshape(   (total_num_system, 1), order='F')
        SP_tem_tmp    = data['SP_tem'].reshape(   (total_num_system, 1), order='F')
        SP_eta_tmp    = data['SP_eta'].reshape(   (total_num_system, 1), order='F')
        SP_sm_rho_tmp = data['SP_sm_rho'].reshape((total_num_system, 1), order='F')
        SP_rho_G_tmp  = data['SP_rho_G'].reshape( (total_num_system, 1), order='F')
        par_con       = data['par_con'].reshape(  (total_num_system, 1), order='F')

        SP_read_tmp = np.hstack((SP_xy_tmp, SP_uv_tmp, SP_rho_S_tmp, SP_pre_tmp, \
                                 SP_tem_tmp, SP_eta_tmp, SP_sm_rho_tmp, SP_rho_G_tmp, par_con))
        SP_read[:, :, i] = SP_read_tmp[:, :]

# particle ID
f = open('./../output_setting/initial_state_SP_kind_{}.dat'.format(file_name), 'rb')
data_type = np.dtype([('SP_kind', endian+str(total_num_system)+'i')])
data = np.fromfile(f, dtype=data_type)
SP_kind = data['SP_kind'].reshape((total_num_system, 1), order='F')

PS_scheme_read = np.zeros((total_num_system, 2, num_save_file+1))
for i in range(num_save_file+1):
        file_num = str(i)
        f = open('./../save_file/PS_scheme_save/PS_scheme_save_{}_{}.dat'.format(file_name, file_num), 'rb')
        data_type = np.dtype([('PS_shift_xy', endian+str(2*total_num_system)+'d')])
        data = np.fromfile(f, dtype=data_type)
        PS_shift_xy_tmp  = data['PS_shift_xy'].reshape((total_num_system, 2), order='F')
        PS_scheme_read[:, :, i] = PS_shift_xy_tmp

print('[Message: 1/8] Reading data has done.')

# translate axes
tra_dis = wall_thickness * delta_x
SP_read[:, 0, :] = SP_read[:, 0, :] - tra_dis
SP_read[:, 1, :] = SP_read[:, 1, :] - tra_dis

# categorize
SP_inner = np.zeros((total_num, 11, num_save_file+1))
SP_wall  = np.zeros((num_wall,  11, num_save_file+1))
PS_shift_inner = np.zeros((total_num, 2, num_save_file+1))
for i in range(num_save_file+1):
        count_inner = 0
        count_wall  = 0
        for j in range(total_num_system):
                if (SP_kind[j] == 0):
                        SP_inner[count_inner, :, i] = SP_read[j, :, i]
                        PS_shift_inner[count_inner, :, i] = PS_scheme_read[j, :, i]
                        count_inner += 1
                else:
                        SP_wall[count_wall, :, i] = SP_read[j, :, i]
                        count_wall += 1

# pass_time
pass_time = np.loadtxt('./../save_file/time_save_{}.dat'.format(file_name))
# # time_step
# time_step = np.loadtxt('./../save_file/time_step_{}.dat'.format(file_name))

print('[Message: 2/8] Sorting data has done.')
read_time = time.perf_counter()
print('[Message     ] Reading & sorting time : {:.2f} [s]'.format(read_time-start_time))
#---------------------
# plot_without_wall
#---------------------
if (plot_without_wall == 1):
    for i in range(num_save_file+1):
        # effective step
        eff_step = i * write_step

        # figure and axis environment
        fig, axs = plt.subplots(1, 1, figsize=(6, 7), facecolor='white', subplot_kw={'facecolor':'white'})
        # margin between figures
        plt.subplots_adjust(left=0.2, right=0.86, bottom=0.02, top=0.9, wspace=0.4, hspace=0.1)

        # plotting
        if (plot_phy == 2):
            fig_1 = axs.scatter(SP_inner[:,0,i]/length, SP_inner[:,1,i]/length, s=SP_size, c=SP_inner[:,2,i]/V_nor, cmap='jet', vmin=-V_max_plot, vmax=V_max_plot)
            bar_1 = plt.colorbar(fig_1, aspect=60, ax=axs, orientation='horizontal', pad=0.18)
            bar_1.set_label(r'$u^{*}$', size=18, labelpad=12) # colorbar label
            bar_1.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)
        if (plot_phy == 3):
            fig_1 = axs.scatter(SP_inner[:,0,i]/length, SP_inner[:,1,i]/length, s=SP_size, c=SP_inner[:,3,i]/V_nor, cmap='jet', vmin=-V_max_plot, vmax=V_max_plot)
            bar_1 = plt.colorbar(fig_1, aspect=60, ax=axs, orientation='horizontal', pad=0.18)
            bar_1.set_label(r'$v^{*}$', size=18, labelpad=12) # colorbar label
            bar_1.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)
        if (plot_phy == 4):
            fig_1 = axs.scatter(SP_inner[:,0,i]/length, SP_inner[:,1,i]/length, s=SP_size, c=SP_inner[:,4,i]/rho_0, cmap='jet', vmin=0.9, vmax=1.1)
            bar_1 = plt.colorbar(fig_1, aspect=60, ax=axs, orientation='horizontal', pad=0.18)
            bar_1.set_label(r'$\rho^{S*}$', size=18, labelpad=12) # colorbar label
            bar_1.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)
        if (plot_phy == 5):
            fig_1 = axs.scatter(SP_inner[:,0,i]/length, SP_inner[:,1,i]/length, s=SP_size, c=SP_inner[:,5,i]/pre_nor, cmap='jet')
            bar_1 = plt.colorbar(fig_1, aspect=60, ax=axs, orientation='horizontal', pad=0.18)
            bar_1.set_label(r'$p^{*}$', size=18, labelpad=12) # colorbar label
            bar_1.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)   
        if (plot_phy == 6):
            fig_1 = axs.scatter(SP_inner[:,0,i]/length, SP_inner[:,1,i]/length, s=SP_size, c=(SP_inner[:,6,i]-T_0)/delta_tem, cmap='jet', vmin=-0.501, vmax=0.501)
            bar_1 = plt.colorbar(fig_1, aspect=60, ax=axs, orientation='horizontal', pad=0.18)
            bar_1.set_label(r'$T^{*}$', size=18, labelpad=12) # colorbar label
            bar_1.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)
        if (plot_phy == 7):
            fig_1 = axs.scatter(SP_inner[:,0,i]/length, SP_inner[:,1,i]/length, s=SP_size, c=SP_inner[:,7,i]/eta_0, cmap='jet')
            bar_1 = plt.colorbar(fig_1, aspect=60, ax=axs, orientation='horizontal', pad=0.18)
            bar_1.set_label(r'$\eta^{*}$', size=18, labelpad=12) # colorbar label
            bar_1.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)
        if (plot_phy == 8):
            fig_1 = axs.scatter(SP_inner[:,0,i]/length, SP_inner[:,1,i]/length, s=SP_size, c=SP_inner[:,8,i]/rho_0, cmap='jet', vmin=0.9, vmax=1.1)
            bar_1 = plt.colorbar(fig_1, aspect=60, ax=axs, orientation='horizontal', pad=0.18)
            bar_1.set_label(r'$\bar{\rho}^{S*}$', size=18, labelpad=12) # colorbar label
            bar_1.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16) 
        if (plot_phy == 9):
            fig_1 = axs.scatter(SP_inner[:,0,i]/length, SP_inner[:,1,i]/length, s=SP_size, c=(SP_inner[:,9,i] + rho_0)/rho_0, cmap='jet')
            bar_1 = plt.colorbar(fig_1, aspect=60, ax=axs, orientation='horizontal', pad=0.18)
            bar_1.set_label(r'$\rho^{G*}$', size=18, labelpad=12) # colorbar label
            bar_1.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)  
        
        # axis labels
        axs.set_xlabel(r'$x^{*}$', fontsize=24, labelpad=14)
        axs.set_ylabel(r'$y^{*}$', fontsize=24, labelpad=20)
        # title
        axs.set_title('{:12.3e}s  ({} step)'.format(pass_time[i], eff_step), fontsize=18, color='k', y=1.04)
        # direction and width of ticks
        axs.tick_params(axis='both', which='major', direction='out', length=3, width=0.8, labelsize=16)
        # width of outer frame
        axs.spines["bottom"].set_linewidth(1.2)
        axs.spines["top"].set_linewidth(1.2)
        axs.spines["right"].set_linewidth(1.2)
        axs.spines["left"].set_linewidth(1.2)
        # lines
        axs.axvline(x=0, ymin=0, ymax=1, color='silver', linewidth=0.8, linestyle='--')
        axs.axvline(x=1, ymin=0, ymax=1, color='silver', linewidth=0.8, linestyle='--')
        axs.axhline(y=0, xmin=0, xmax=1, color='silver', linewidth=0.8, linestyle='--')
        axs.axhline(y=1, xmin=0, xmax=1, color='silver', linewidth=0.8, linestyle='--')
        # lim
        axs.set_xlim(0-0.05, 1+0.05)
        axs.set_ylim(0-0.05, 1+0.05)
        # save       
        fig.savefig('./../fig/movie_material/SP_{}_{}_{}.png'.format(file_name, plot_phy, i), format='png', dpi=300, transparent=False)
        # close
        plt.close()

    # print
    print('[Message: 3/8] Plot_movie_material     without wall has     done.')
else:
    print('[Message: 3/8] Plot_movie_material     without wall has NOT done.')
#---------------------
# plot_without_wall
#---------------------
if (plot_with_wall == 1):
    for i in range(num_save_file+1):
        # effective step
        eff_step = i * write_step

        # figure and axis environment
        fig, axs = plt.subplots(1, 1, figsize=(6, 7), facecolor='white', subplot_kw={'facecolor':'white'})
        # margin between figures
        plt.subplots_adjust(left=0.2, right=0.86, bottom=0.02, top=0.9, wspace=0.4, hspace=0.1)

        # plotting
        if (plot_phy == 2):
            fig_1 = axs.scatter(SP_read[:,0,i]/length, SP_read[:,1,i]/length, s=SP_size, c=SP_read[:,2,i]/V_nor, cmap='jet', vmin=-V_max_plot, vmax=V_max_plot)
            bar_1 = plt.colorbar(fig_1, aspect=60, ax=axs, orientation='horizontal', pad=0.18)
            bar_1.set_label(r'$u^{*}$', size=18, labelpad=12) # colorbar label
            bar_1.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)
        if (plot_phy == 3):
            fig_1 = axs.scatter(SP_read[:,0,i]/length, SP_read[:,1,i]/length, s=SP_size, c=SP_read[:,3,i]/V_nor, cmap='jet', vmin=-V_max_plot, vmax=V_max_plot)
            bar_1 = plt.colorbar(fig_1, aspect=60, ax=axs, orientation='horizontal', pad=0.18)
            bar_1.set_label(r'$v^{*}$', size=18, labelpad=12) # colorbar label
            bar_1.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)
        if (plot_phy == 4):
            fig_1 = axs.scatter(SP_read[:,0,i]/length, SP_read[:,1,i]/length, s=SP_size, c=SP_read[:,4,i]/rho_0, cmap='jet', vmin=0.9, vmax=1.1)
            bar_1 = plt.colorbar(fig_1, aspect=60, ax=axs, orientation='horizontal', pad=0.18)
            bar_1.set_label(r'$\rho^{S*}$', size=18, labelpad=12) # colorbar label
            bar_1.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)
        if (plot_phy == 5):
            fig_1 = axs.scatter(SP_read[:,0,i]/length, SP_read[:,1,i]/length, s=SP_size, c=SP_read[:,5,i]/pre_nor, cmap='jet')
            bar_1 = plt.colorbar(fig_1, aspect=60, ax=axs, orientation='horizontal', pad=0.18)
            bar_1.set_label(r'$p^{*}$', size=18, labelpad=12) # colorbar label
            bar_1.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)   
        if (plot_phy == 6):
            fig_1 = axs.scatter(SP_read[:,0,i]/length, SP_read[:,1,i]/length, s=SP_size, c=(SP_read[:,6,i]-T_0)/delta_tem, cmap='jet', vmin=-0.501, vmax=0.501)
            bar_1 = plt.colorbar(fig_1, aspect=60, ax=axs, orientation='horizontal', pad=0.18)
            bar_1.set_label(r'$T^{*}$', size=18, labelpad=12) # colorbar label
            bar_1.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)
        if (plot_phy == 7):
            fig_1 = axs.scatter(SP_read[:,0,i]/length, SP_read[:,1,i]/length, s=SP_size, c=SP_read[:,7,i]/eta_0, cmap='jet')
            bar_1 = plt.colorbar(fig_1, aspect=60, ax=axs, orientation='horizontal', pad=0.18)
            bar_1.set_label(r'$\eta^{*}$', size=18, labelpad=12) # colorbar label
            bar_1.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)
        if (plot_phy == 8):
            fig_1 = axs.scatter(SP_read[:,0,i]/length, SP_read[:,1,i]/length, s=SP_size, c=SP_read[:,8,i]/rho_0, cmap='jet', vmin=0.9, vmax=1.1)
            bar_1 = plt.colorbar(fig_1, aspect=60, ax=axs, orientation='horizontal', pad=0.18)
            bar_1.set_label(r'$\bar{\rho}^{S*}$', size=18, labelpad=12) # colorbar label
            bar_1.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16) 
        if (plot_phy == 9):
            fig_1 = axs.scatter(SP_read[:,0,i]/length, SP_read[:,1,i]/length, s=SP_size, c=(SP_read[:,9,i] + rho_0)/rho_0, cmap='jet')
            bar_1 = plt.colorbar(fig_1, aspect=60, ax=axs, orientation='horizontal', pad=0.18)
            bar_1.set_label(r'$\rho^{G*}$', size=18, labelpad=12) # colorbar label
            bar_1.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)  
        
        # axis labels
        axs.set_xlabel(r'$x^{*}$', fontsize=24, labelpad=14)
        axs.set_ylabel(r'$y^{*}$', fontsize=24, labelpad=20)
        # title
        axs.set_title('{:12.3e}s  ({} step)'.format(pass_time[i], eff_step), fontsize=18, color='k', y=1.04)
        # direction and width of ticks
        axs.tick_params(axis='both', which='major', direction='out', length=3, width=0.8, labelsize=16)
        # width of outer frame
        axs.spines["bottom"].set_linewidth(1.2)
        axs.spines["top"].set_linewidth(1.2)
        axs.spines["right"].set_linewidth(1.2)
        axs.spines["left"].set_linewidth(1.2)
        # lines
        axs.axvline(x=0, ymin=0, ymax=1, color='silver', linewidth=0.8, linestyle='--')
        axs.axvline(x=1, ymin=0, ymax=1, color='silver', linewidth=0.8, linestyle='--')
        axs.axhline(y=0, xmin=0, xmax=1, color='silver', linewidth=0.8, linestyle='--')
        axs.axhline(y=1, xmin=0, xmax=1, color='silver', linewidth=0.8, linestyle='--')
        # lim
        axs.set_xlim(0-0.18, 1+0.18)
        axs.set_ylim(0-0.18, 1+0.18)
        # save       
        fig.savefig('./../fig/movie_material_wall/SP_wall_{}_{}_{}.png'.format(file_name, plot_phy, i), format='png', dpi=300, transparent=False)
        # close
        plt.close()

    # print
    print('[Message: 4/8] Plot_movie_material     with    wall has     done.')
else:
    print('[Message: 4/8] Plot_movie_material     with    wall has NOT done.')

#----------------------
# plot_all_witout_wall
#----------------------
if (plot_all_without_wall == 1):
    for i in range(num_save_file+1):
        # effective step
        eff_step = i * write_step
        
        # figure and axis environment
        fig, axs = plt.subplots(3, 3, figsize=(18, 22), facecolor='white', subplot_kw={'facecolor':'white'})
        
        # margin between figure
        plt.subplots_adjust(left=0.08, right=0.95, bottom=0.02, top=0.95, wspace=0.35, hspace=0.1)
        
        # u
        fig_1 = axs[0, 0].scatter(SP_inner[:,0,i]/length, SP_inner[:,1,i]/length, s=SP_size, c=SP_inner[:,2,i]/V_nor, cmap='jet', vmin=-V_max_plot, vmax=V_max_plot)
        bar_1 = plt.colorbar(fig_1, aspect=60, ax=axs[0, 0], orientation='horizontal', pad=0.18)
        bar_1.set_label(r'$u^{*}$', size=24, labelpad=14)
        bar_1.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=18)

        # v
        fig_2 = axs[0, 1].scatter(SP_inner[:,0,i]/length, SP_inner[:,1,i]/length, s=SP_size, c=SP_inner[:,3,i]/V_nor, cmap='jet', vmin=-V_max_plot, vmax=V_max_plot)
        bar_2 = plt.colorbar(fig_2, aspect=60, ax=axs[0, 1], orientation='horizontal', pad=0.18)
        bar_2.set_label(r'$v^{*}$', size=24, labelpad=14)
        bar_2.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=18) 

        # SPH density
        fig_3 = axs[0, 2].scatter(SP_inner[:,0,i]/length, SP_inner[:,1,i]/length, s=SP_size, c=SP_inner[:,4,i]/rho_0, cmap='jet', vmin=0.9, vmax=1.1)
        bar_3 = plt.colorbar(fig_3, aspect=60, ax=axs[0, 2], orientation='horizontal', pad=0.18)
        bar_3.set_label(r'$\rho^{S*}$', size=24, labelpad=14)
        bar_3.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=18) 

        # pressure
        fig_4 = axs[1, 0].scatter(SP_inner[:,0,i]/length, SP_inner[:,1,i]/length, s=SP_size, c=SP_inner[:,5,i]/pre_nor, cmap='jet')
        bar_4 = plt.colorbar(fig_4, aspect=60, ax=axs[1, 0], orientation='horizontal', pad=0.18)
        bar_4.set_label(r'$p^{*}$', size=24, labelpad=14)
        bar_4.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=18)

        # temperature
        fig_5 = axs[1, 1].scatter(SP_inner[:,0,i]/length, SP_inner[:,1,i]/length, s=SP_size, c=(SP_inner[:,6,i]-T_0)/delta_tem, cmap='jet', vmin=-0.501, vmax=0.501)
        bar_5 = plt.colorbar(fig_5, aspect=60, ax=axs[1, 1], orientation='horizontal', pad=0.18)
        bar_5.set_label(r'$T^{*}$', size=24, labelpad=14)
        bar_5.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=18)

        # eta
        fig_6 = axs[1, 2].scatter(SP_inner[:,0,i]/length, SP_inner[:,1,i]/length, s=SP_size, c=SP_inner[:,7,i]/eta_0, cmap='jet')
        bar_6 = plt.colorbar(fig_6, aspect=60, ax=axs[1, 2], orientation='horizontal', pad=0.18)
        bar_6.set_label(r'$\eta^{*}$', size=24, labelpad=14)
        bar_6.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=18)

        # smoothed density
        fig_7 = axs[2, 0].scatter(SP_inner[:,0,i]/length, SP_inner[:,1,i]/length, s=SP_size, c=SP_inner[:,8,i]/rho_0, cmap='jet', vmin=0.9, vmax=1.1)
        bar_7 = plt.colorbar(fig_7, aspect=60, ax=axs[2, 0], orientation='horizontal', pad=0.18)
        bar_7.set_label(r'$\bar{\rho}^{S*}$', size=24, labelpad=14)
        bar_7.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=18)

        # gravitational density
        fig_8 = axs[2, 1].scatter(SP_inner[:,0,i]/length, SP_inner[:,1,i]/length, s=SP_size, c=(SP_inner[:,9,i] + rho_0)/rho_0, cmap='jet')
        bar_8 = plt.colorbar(fig_8, aspect = 60, ax = axs[2, 1], orientation = 'horizontal', pad = 0.18)
        bar_8.set_label(r'$\rho^{G*}$', size=24, labelpad=14)
        bar_8.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=18)

        # particle consistency
        fig_9 = axs[2, 2].scatter(SP_inner[:,0,i]/length, SP_inner[:,1,i]/length, s=SP_size, c=SP_inner[:,10,i], cmap='jet', vmin=0.9, vmax=1.1)
        bar_9 = plt.colorbar(fig_9, aspect=60, ax=axs[2, 2], orientation='horizontal', pad=0.18)
        bar_9.set_label(r'$\varpi$', size=24, labelpad=14)
        bar_9.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=18)

        # axis labels, ticks, outer frame, lines, lim
        for j in range(3):
            for k in range(3):
                    axs[j, k].set_xlabel(r'$x^{*}$', fontsize=24, labelpad=18)
                    axs[j, k].set_ylabel(r'$y^{*}$', fontsize=24, labelpad=18)
                    axs[j, k].tick_params(axis='both', which='major', direction='out', length=3, width=0.8, labelsize=18)
                    axs[j, k].spines["bottom"].set_linewidth(1.2)
                    axs[j, k].spines["top"].set_linewidth(1.2)
                    axs[j, k].spines["right"].set_linewidth(1.2)
                    axs[j, k].spines["left"].set_linewidth(1.2)
                    axs[j, k].axvline(x=0, ymin=0, ymax=1, color='silver', linewidth=0.8, linestyle='--')
                    axs[j, k].axvline(x=1, ymin=0, ymax=1, color='silver', linewidth=0.8, linestyle='--')
                    axs[j, k].axhline(y=0, xmin=0, xmax=1, color='silver', linewidth=0.8, linestyle='--')
                    axs[j, k].axhline(y=1, xmin=0, xmax=1, color='silver', linewidth=0.8, linestyle='--')
                    axs[j, k].set_xlim(0-0.05, 1+0.05)
                    axs[j, k].set_ylim(0-0.05, 1+0.05)
        # title
        plt.suptitle('{:12.3e}s  ({}step)'.format(pass_time[i], eff_step), fontsize=30, color='k', x=0.515, y=0.985)
        # save
        fig.savefig('./../fig/movie_material_all/SP_all_{}_{}.png'.format(file_name, i), format='png', dpi=300, transparent=False)
        # close
        plt.close()

    print('[Message: 5/8] Plot_movie_meterial_all without wall has     done.')
else:
    print('[Message: 5/8] Plot_movie_meterial_all without wall has NOT done.')

#----------------------
# plot_all_with_wall
#----------------------
if (plot_all_with_wall == 1):
    for i in range(num_save_file+1):
        # effective step
        eff_step = i * write_step
        
        # figure and axis environment
        fig, axs = plt.subplots(3, 3, figsize=(18, 22), facecolor='white', subplot_kw={'facecolor':'white'})
        
        # margin between figure
        plt.subplots_adjust(left=0.08, right=0.95, bottom=0.02, top=0.95, wspace=0.35, hspace=0.1)
        
        # u
        fig_1 = axs[0, 0].scatter(SP_read[:,0,i]/length, SP_read[:,1,i]/length, s=SP_size, c=SP_read[:,2,i]/V_nor, cmap='jet', vmin=-V_max_plot, vmax=V_max_plot)
        bar_1 = plt.colorbar(fig_1, aspect=60, ax=axs[0, 0], orientation='horizontal', pad=0.18)
        bar_1.set_label(r'$u^{*}$', size=24, labelpad=14)
        bar_1.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=18)

        # v
        fig_2 = axs[0, 1].scatter(SP_read[:,0,i]/length, SP_read[:,1,i]/length, s=SP_size, c=SP_read[:,3,i]/V_nor, cmap='jet', vmin=-V_max_plot, vmax=V_max_plot)
        bar_2 = plt.colorbar(fig_2, aspect=60, ax=axs[0, 1], orientation='horizontal', pad=0.18)
        bar_2.set_label(r'$v^{*}$', size=24, labelpad=14)
        bar_2.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=18) 

        # SPH density
        fig_3 = axs[0, 2].scatter(SP_read[:,0,i]/length, SP_read[:,1,i]/length, s=SP_size, c=SP_read[:,4,i]/rho_0, cmap='jet', vmin=0.9, vmax=1.1)
        bar_3 = plt.colorbar(fig_3, aspect=60, ax=axs[0, 2], orientation='horizontal', pad=0.18)
        bar_3.set_label(r'$\rho^{S*}$', size=24, labelpad=14)
        bar_3.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=18) 

        # pressure
        fig_4 = axs[1, 0].scatter(SP_read[:,0,i]/length, SP_read[:,1,i]/length, s=SP_size, c=SP_read[:,5,i]/pre_nor, cmap='jet')
        bar_4 = plt.colorbar(fig_4, aspect=60, ax=axs[1, 0], orientation='horizontal', pad=0.18)
        bar_4.set_label(r'$p^{*}$', size=24, labelpad=14)
        bar_4.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=18)

        # temperature
        fig_5 = axs[1, 1].scatter(SP_read[:,0,i]/length, SP_read[:,1,i]/length, s=SP_size, c=(SP_read[:,6,i]-T_0)/delta_tem, cmap='jet', vmin=-0.501, vmax=0.501)
        bar_5 = plt.colorbar(fig_5, aspect=60, ax=axs[1, 1], orientation='horizontal', pad=0.18)
        bar_5.set_label(r'$T^{*}$', size=24, labelpad=14)
        bar_5.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=18)

        # eta
        fig_6 = axs[1, 2].scatter(SP_read[:,0,i]/length, SP_read[:,1,i]/length, s=SP_size, c=SP_read[:,7,i]/eta_0, cmap='jet')
        bar_6 = plt.colorbar(fig_6, aspect=60, ax=axs[1, 2], orientation='horizontal', pad=0.18)
        bar_6.set_label(r'$\eta^{*}$', size=24, labelpad=14)
        bar_6.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=18)

        # smoothed density
        fig_7 = axs[2, 0].scatter(SP_read[:,0,i]/length, SP_read[:,1,i]/length, s=SP_size, c=SP_read[:,8,i]/rho_0, cmap='jet', vmin=0.9, vmax=1.1)
        bar_7 = plt.colorbar(fig_7, aspect=60, ax=axs[2, 0], orientation='horizontal', pad=0.18)
        bar_7.set_label(r'$\bar{\rho}^{S*}$', size=24, labelpad=14)
        bar_7.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=18)

        # gravitational density
        fig_8 = axs[2, 1].scatter(SP_read[:,0,i]/length, SP_read[:,1,i]/length, s=SP_size, c=(SP_read[:,9,i] + rho_0)/rho_0, cmap='jet')
        bar_8 = plt.colorbar(fig_8, aspect = 60, ax = axs[2, 1], orientation = 'horizontal', pad = 0.18)
        bar_8.set_label(r'$\rho^{G*}$', size=24, labelpad=14)
        bar_8.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=18)

        # particle consistency
        fig_9 = axs[2, 2].scatter(SP_read[:,0,i]/length, SP_read[:,1,i]/length, s=SP_size, c=SP_read[:,10,i], cmap='jet', vmin=0.9, vmax=1.1)
        bar_9 = plt.colorbar(fig_9, aspect=60, ax=axs[2, 2], orientation='horizontal', pad=0.18)
        bar_9.set_label(r'$\varpi$', size=24, labelpad=14)
        bar_9.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=18)

        # axis labels, ticks, outer frame, lines
        for j in range(3):
            for k in range(3):
                    axs[j, k].set_xlabel(r'$x^{*}$', fontsize=24, labelpad=18)
                    axs[j, k].set_ylabel(r'$y^{*}$', fontsize=24, labelpad=18)
                    axs[j, k].tick_params(axis='both', which='major', direction='out', length=3, width=0.8, labelsize=18)
                    axs[j, k].spines["bottom"].set_linewidth(1.2)
                    axs[j, k].spines["top"].set_linewidth(1.2)
                    axs[j, k].spines["right"].set_linewidth(1.2)
                    axs[j, k].spines["left"].set_linewidth(1.2)
                    axs[j, k].axvline(x=0, ymin=0, ymax=1, color='silver', linewidth=0.8, linestyle='--')
                    axs[j, k].axvline(x=1, ymin=0, ymax=1, color='silver', linewidth=0.8, linestyle='--')
                    axs[j, k].axhline(y=0, xmin=0, xmax=1, color='silver', linewidth=0.8, linestyle='--')
                    axs[j, k].axhline(y=1, xmin=0, xmax=1, color='silver', linewidth=0.8, linestyle='--')
                    axs[j, k].set_xlim(0-0.18, 1+0.18)
                    axs[j, k].set_ylim(0-0.18, 1+0.18)
        # title
        plt.suptitle('{:12.3e}s  ({}step)'.format(pass_time[i], eff_step), fontsize=30, color='k', x=0.515, y=0.985)
        # save
        fig.savefig('./../fig/movie_material_all_wall/SP_all_wall_{}_{}.png'.format(file_name, i), format='png', dpi=300, transparent=False)
        # close
        plt.close()

    print('[Message: 6/8] Plot_movie_meterial_all with    wall has     done.')
else:
    print('[Message: 6/8] Plot_movie_meterial_all with    wall has NOT done.')

#---------------------
# PS_shift
#---------------------
if (plot_PS_shift == 1):
    for i in range(1, num_save_file+1):
        # effective step
        eff_step = i * write_step

        # figure and axis environment
        fig, axs = plt.subplots(1, 3, figsize=(18, 8), facecolor='white', subplot_kw={'facecolor':'white'})
                
        # margin between figure
        plt.subplots_adjust(left=0.08, right=0.95, bottom=0.02, top=0.88, wspace=0.3, hspace=0.1)
                
        # shift x
        fig_1 = axs[0].scatter(SP_inner[:,0,i]/length, SP_inner[:,1,i]/length, s=SP_size, \
                c=PS_shift_inner[:,0,i]/np.abs((SP_inner[:,2,i]*time_step)), cmap='jet', vmin=-0.2, vmax=0.2)
        bar_1 = plt.colorbar(fig_1, aspect=60, ax=axs[0], orientation='horizontal', pad=0.18)
        bar_1.set_label(r'$\delta x^{*}$', size=28, labelpad=14)
        bar_1.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=24)
        axs[0].quiver(SP_inner[:,0,i]/length, SP_inner[:,1,i]/length, PS_shift_inner[:,0,i], np.zeros((total_num)), color='k')

        # shift y
        fig_2 = axs[1].scatter(SP_inner[:,0,i]/length, SP_inner[:,1,i]/length, s=SP_size, \
                c=PS_shift_inner[:,1,i]/np.abs((SP_inner[:,3,i]*time_step)), cmap='jet', vmin=-0.2, vmax=0.2)
        bar_2 = plt.colorbar(fig_2, aspect=60, ax=axs[1], orientation='horizontal', pad=0.18)
        bar_2.set_label(r'$\delta y^{*}$', size=28, labelpad=14)
        bar_2.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=24)
        axs[1].quiver(SP_inner[:,0,i]/length, SP_inner[:,1,i]/length, np.zeros((total_num)), PS_shift_inner[:,1,i], color='k')

        # shift dis
        fig_3 = axs[2].scatter(SP_inner[:,0,i]/length, SP_inner[:,1,i]/length, s=SP_size, \
                c=np.sqrt(PS_shift_inner[:,0,i]**2+PS_shift_inner[:,1,i]**2) \
                /(np.sqrt((SP_inner[:,2,i]**2+SP_inner[:,3,i]**2))*time_step), cmap='jet', vmin=0.0, vmax=0.2)
        bar_3 = plt.colorbar(fig_3, aspect=60, ax=axs[2], orientation='horizontal', pad=0.18)
        bar_3.set_label(r'$\delta r^{*}$', size=28, labelpad=14)
        bar_3.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=24) 
        axs[2].quiver(SP_inner[:,0,i]/length, SP_inner[:,1,i]/length, PS_shift_inner[:,0,i], PS_shift_inner[:,1,i], color='k')

        axs[0].set_ylabel(r'$y^{*}$', fontsize=30, labelpad=20)
        # axis labels, ticks, outer frame, lines, lim
        for j in range(3):
                axs[j].set_xlabel(r'$x^{*}$', fontsize=30, labelpad=14)
                axs[j].tick_params(axis='both', which='major', direction='out', length=3, width=0.8, labelsize=24)
                axs[j].spines["bottom"].set_linewidth(1.2)
                axs[j].spines["top"].set_linewidth(1.2)
                axs[j].spines["right"].set_linewidth(1.2)
                axs[j].spines["left"].set_linewidth(1.2)
                axs[j].axvline(x=0, ymin=0, ymax=1, color='silver', linewidth=0.8, linestyle='--')
                axs[j].axvline(x=1, ymin=0, ymax=1, color='silver', linewidth=0.8, linestyle='--')
                axs[j].axhline(y=0, xmin=0, xmax=1, color='silver', linewidth=0.8, linestyle='--')
                axs[j].axhline(y=1, xmin=0, xmax=1, color='silver', linewidth=0.8, linestyle='--')
                axs[j].set_xlim(0-0.05, 1+0.05)
                axs[j].set_ylim(0-0.05, 1+0.05)
        # title
        plt.suptitle('{:12.3e}s  ({}step)'.format(pass_time[i], eff_step), fontsize=30, color='k', x=0.515, y=0.975)
        # save
        fig.savefig('./../fig/PS_shift/PS_shift_{}_{}.png'.format(file_name, i), format='png', dpi=300, transparent=False)
        # close
        plt.close()
    print('[Message: 7/8] Plot_PS_shift                        has     done.')
else:
    print('[Message: 7/8] Plot_PS_shift                        has NOT done.')


print('[Message: 8/8] Plot_movie_material                      has completed !')
end_time = time.perf_counter()
print('-----------------------------------------------------')
print('[Message] Total time : {:.2f} [s]'.format(end_time-start_time))
print('-----------------------------------------------------')
print('')