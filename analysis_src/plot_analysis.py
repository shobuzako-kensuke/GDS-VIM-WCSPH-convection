#==========================================================================================#
#                           thermal convection using GDS-WCSPH                             #
#------------------------------------------------------------------------------------------#
#                          Copyright by Kensuke Shobuzako (2022)                           #
#==========================================================================================#

#=====================
# CHARM
#=====================
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
print('[Message] plot_analysis has started.')
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
                              ('SP_rho_S',  endian+str(  total_num_system)+'d'), \
                              ('SP_pre',    endian+str(  total_num_system)+'d'), \
                              ('SP_tem',    endian+str(  total_num_system)+'d'), \
                              ('SP_eta',    endian+str(  total_num_system)+'d'), \
                              ('SP_sm_rho', endian+str(  total_num_system)+'d'), \
                              ('SP_rho_G',  endian+str(  total_num_system)+'d'), \
                              ('par_con',   endian+str(  total_num_system)+'d')])
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

# EOM
EOM_read = np.zeros((total_num_system, 7, num_save_file+1)) # num_save_file = initial_state + save_state
for i in range(num_save_file+1):
        file_num = str(i)
        f = open('./../save_file/EOM_save/EOM_save_{}_{}.dat'.format(file_name, file_num), 'rb')
        data_type = np.dtype([('vis_EOM',     endian+str(2*total_num_system)+'d'), \
                              ('buo_EOM',     endian+str(  total_num_system)+'d'), \
                              ('pre_EOM',     endian+str(2*total_num_system)+'d'), \
                              ('ine_EOM',     endian+str(2*total_num_system)+'d')])
        data = np.fromfile(f, dtype=data_type)
        vis_EOM_tmp = data['vis_EOM'].reshape((total_num_system, 2), order='F')
        buo_EOM_tmp = data['buo_EOM'].reshape((total_num_system, 1), order='F')
        pre_EOM_tmp = data['pre_EOM'].reshape((total_num_system, 2), order='F')
        ine_EOM_tmp = data['ine_EOM'].reshape((total_num_system, 2), order='F')

        EOM_read_tmp = np.hstack((vis_EOM_tmp, buo_EOM_tmp, pre_EOM_tmp, ine_EOM_tmp))
        EOM_read[:, :, i] = EOM_read_tmp[:, :]

# particle ID
f = open('./../output_setting/initial_state_SP_kind_{}.dat'.format(file_name), 'rb')
data_type = np.dtype([('SP_kind', endian+str(total_num_system)+'i')])
data = np.fromfile(f, dtype=data_type)
SP_kind = data['SP_kind'].reshape((total_num_system, 1), order='F')

# MP
MP_read = np.zeros((num_MP, 7, num_save_file+1))
for i in range(num_save_file+1):
        file_num = str(i)
        f = open('./../save_file/MP_save/MP_save_{}_{}.dat'.format(file_name, file_num), 'rb')
        data_type = np.dtype([('MP_xy',  endian+str(2*num_MP)+'d'), \
                              ('MP_uv',  endian+str(2*num_MP)+'d'), \
                              ('MP_pre', endian+str(num_MP)+'d'), \
                              ('MP_tem', endian+str(num_MP)+'d'), \
                              ('MP_NU',  endian+str(num_MP)+'d') ])
        data = np.fromfile(f, dtype=data_type)
        MP_xy_tmp  = data['MP_xy'].reshape( (num_MP, 2), order='F')
        MP_uv_tmp  = data['MP_uv'].reshape( (num_MP, 2), order='F')
        MP_pre_tmp = data['MP_pre'].reshape((num_MP, 1), order='F')
        MP_tem_tmp = data['MP_tem'].reshape((num_MP, 1), order='F')
        MP_NU_tmp  = data['MP_NU'].reshape( (num_MP, 1), order='F')

        MP_read_tmp = np.hstack((MP_xy_tmp, MP_uv_tmp, MP_pre_tmp, MP_tem_tmp, MP_NU_tmp))
        MP_read[:, :, i] = MP_read_tmp[:, :]

# V_rms_ave
V_rms_ave_read = np.loadtxt('./../save_file/V_rms_save_{}.dat'.format(file_name))
V_rms_ave = np.zeros((num_save_file+1))
# pass_time
pass_time_read = np.loadtxt('./../save_file/time_save_{}.dat'.format(file_name))
pass_time = np.zeros((num_save_file+1))
# # time_step
# time_step_read = np.loadtxt('./../save_file/time_step_{}.dat'.format(file_name))
# time_step = np.zeros((num_save_file+1))
# # RSS_arg
# RSS_arg_read = np.loadtxt('./../save_file/RSS_{}.dat'.format(file_name))
# RSS_arg = np.zeros((num_save_file+1))

for i in range(num_save_file+1):
        V_rms_ave[i] = V_rms_ave_read[i]
        pass_time[i] = pass_time_read[i]
        # time_step[i] = time_step_read[i]
        # RSS_arg[i] = RSS_arg_read[i]

# PS_scheme
PS_scheme_read = np.zeros((total_num_system, 2, num_save_file+1))
for i in range(num_save_file+1):
        file_num = str(i)
        f = open('./../save_file/PS_scheme_save/PS_scheme_save_{}_{}.dat'.format(file_name, file_num), 'rb')
        data_type = np.dtype([('PS_shift_xy', endian+str(2*total_num_system)+'d')])
        data = np.fromfile(f, dtype=data_type)
        PS_shift_xy_tmp  = data['PS_shift_xy'].reshape((total_num_system, 2), order='F')
        PS_scheme_read[:, :, i] = PS_shift_xy_tmp

print('[Message:  1/16] Reading data has done.')

# translate axes
tra_dis = wall_thickness * delta_x
SP_read[:, 0, :] = SP_read[:, 0, :] - tra_dis
SP_read[:, 1, :] = SP_read[:, 1, :] - tra_dis
MP_read[:, 0, :] = MP_read[:, 0, :] - tra_dis
MP_read[:, 1, :] = MP_read[:, 1, :] - tra_dis

# categorize
SP_inner  = np.zeros((total_num, 11, num_save_file+1))
EOM_inner = np.zeros((total_num, 7,  num_save_file+1))
SP_wall   = np.zeros((num_wall,  11, num_save_file+1))
MP_b  = np.zeros((num_MP_side, 7, num_save_file+1))
MP_hm = np.zeros((num_MP_side, 7, num_save_file+1))
MP_t  = np.zeros((num_MP_side, 7, num_save_file+1))
MP_l  = np.zeros((num_MP_side, 7, num_save_file+1))
MP_vm = np.zeros((num_MP_side, 7, num_save_file+1))
MP_r  = np.zeros((num_MP_side, 7, num_save_file+1))
PS_shift_inner = np.zeros((total_num, 2, num_save_file+1))

for i in range(num_save_file+1):
        count_inner = 0
        count_wall  = 0
        for j in range(total_num_system):
                if (SP_kind[j] == 0):
                        SP_inner[count_inner, :, i] = SP_read[j, :, i]
                        EOM_inner[count_inner, :, i] = EOM_read[j, :, i]
                        PS_shift_inner[count_inner, :, i] = PS_scheme_read[j, :, i]
                        count_inner += 1
                else:
                        SP_wall[count_wall, :, i] = SP_read[j, :, i]
                        count_wall += 1

for i in range(num_save_file+1):
        for j in range(0*num_MP_side, 1*num_MP_side):
                k = j
                MP_b[k, :, i] = MP_read[j, :, i]
        for j in range(1*num_MP_side, 2*num_MP_side):
                k = j - 1*num_MP_side
                MP_hm[k, :, i] = MP_read[j, :, i]
        for j in range(2*num_MP_side, 3*num_MP_side):
                k = j - 2*num_MP_side
                MP_t[k, :, i] = MP_read[j, :, i]
        for j in range(3*num_MP_side, 4*num_MP_side):
                k = j - 3*num_MP_side
                MP_l[k, :, i] = MP_read[j, :, i]
        for j in range(4*num_MP_side, 5*num_MP_side):
                k = j - 4*num_MP_side
                MP_vm[k, :, i] = MP_read[j, :, i]
        for j in range(5*num_MP_side, 6*num_MP_side):
                k = j - 5*num_MP_side
                MP_r[k, :, i] = MP_read[j, :, i]
                                
print('[Message:  2/16] Sorting data has done.')
read_time = time.perf_counter()
print('[Message       ] Reading & sorting time : {:.2f} [s]'.format(read_time-start_time))

#---------------------
# final state
#---------------------
i = num_save_file
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
# title
plt.suptitle('{:12.3e}s  ({}step)'.format(pass_time[-1], eff_step), fontsize=30, color='k', x=0.515, y=0.985)

# save
fig.savefig('./../fig/analysis/final_state_{}.png'.format(file_name), format='png', dpi=300, transparent=False)
fig.savefig('./../fig/analysis/final_state_{}.pdf'.format(file_name), format='pdf', transparent=True)

# close
plt.close()
print('[Message:  3/16] Plot final_state      has done.')

#---------------------
# V_rms & local Mach
#---------------------
V_rms_final = V_rms_ave[-1] / V_nor
V_arg = np.zeros((num_system, num_save_file+1))
for i in range(num_system):
        for j in range(num_save_file+1):
                V_arg[i, j] = np.sqrt(SP_inner[i, 2, j]**2 + SP_inner[i, 3, j]**2)

V_max_arg = np.zeros((num_save_file+1))
for i in range(num_save_file+1):
        V_max_arg[i] = np.max(V_arg[:, i], axis=0)

Mach_final = V_max_arg[-1] / RSS

# figure and axis environment
fig, axs = plt.subplots(1, 2, figsize=(16, 7.3), facecolor='white', subplot_kw={'facecolor':'white'})

# margin between figures
plt.subplots_adjust(left=0.12, right=0.95, bottom=0.15, top=0.9, wspace=0.5, hspace=0.3)

# plotting
axs[0].plot(pass_time, V_rms_ave/V_nor,   color='k', linewidth=2.5, linestyle='-' )
axs[1].plot(pass_time, V_max_arg/RSS, color='k', linewidth=2.5, linestyle='-' )

# axis labels
axs[0].set_xlabel(r'time $(\mathrm{s})$',    fontsize=26, labelpad=18)
axs[0].set_ylabel(r'$V_{\mathrm{rms}}^{*}$', fontsize=26, labelpad=18)
axs[1].set_xlabel(r'time $(\mathrm{s})$',    fontsize=26, labelpad=18)
axs[1].set_ylabel(r'Mach number',            fontsize=26, labelpad=18)

axs[0].set_title(r'$\left(V_{}^{}\right)_{}={:.5f}$'.format('{\mathrm{rms}}', '*', '{\mathrm{fin}}', V_rms_final), \
        fontsize=26, color='k', y=1.03)
axs[1].set_title(r'$M_{}={:.5f}$'.format('{\mathrm{fin}}', Mach_final), \
        fontsize=26, color='k', y=1.03)

# grid, ticks, outer_frame
for i in range(2):
        axs[i].grid(which='major', color='silver', linewidth=1)
        axs[i].tick_params(axis='both', which='major', direction='out', length=10, width=1.2, labelsize=24)
        axs[i].tick_params(axis='both', which='minor', direction='out', length=10, width=1.2, labelsize=24)
        axs[i].spines["bottom"].set_linewidth(1.2)
        axs[i].spines["top"].set_linewidth(1.2)
        axs[i].spines["right"].set_linewidth(1.2)
        axs[i].spines["left"].set_linewidth(1.2)

# save
fig.savefig('./../fig/analysis/V_rms_local_Mach_{}.png'.format(file_name), format='png', dpi=300, transparent=False)
fig.savefig('./../fig/analysis/V_rms_local_Mach_{}.pdf'.format(file_name), format='pdf', transparent=True)

# close
plt.close()

# print
print('[Message:  4/16] Plot V_rms_local_Mach has done.')

#---------------------
# final V profile
#---------------------
SP_max_arg = np.max(SP_inner[:,:,-1], axis=0)
max_u = SP_max_arg[2] / V_nor
max_v = SP_max_arg[3] / V_nor

# figure and axis environment
fig, axs = plt.subplots(1, 2, figsize=(16, 7.3), facecolor='white', subplot_kw={'facecolor':'white'})

# margin between figures
plt.subplots_adjust(left=0.12, right=0.95, bottom=0.15, top=0.9, wspace=0.5, hspace=0.3)

# plotting
axs[0].plot(MP_vm[:,2,-1]/V_nor,  MP_vm[:,1,-1]/length, linestyle='-', marker='o', c='k')
axs[1].plot(MP_hm[:,0,-1]/length, MP_hm[:,3,-1]/V_nor,  linestyle='-', marker='o', c='k')

# axis labels
axs[0].set_xlabel(r'$u^{*}$', fontsize=32, labelpad=18)
axs[0].set_ylabel(r'$y^{*}$', fontsize=32, labelpad=18)
axs[1].set_xlabel(r'$x^{*}$', fontsize=32, labelpad=18)
axs[1].set_ylabel(r'$v^{*}$', fontsize=32, labelpad=18)

# title
axs[0].set_title(r'$u_{}^{}=~~{:.5f}$'.format('{\max}', '*', max_u), fontsize=26, color='k', y=1.03)
axs[1].set_title(r'$v_{}^{}=~~{:.5f}$'.format('{\max}', '*', max_v), fontsize=26, color='k', y=1.03)
    
# grid, ticks, outer_frame
for i in range(2):
        axs[i].grid(which='major', color='silver', linewidth=1)
        axs[i].tick_params(axis='both', which='major', direction='out', length=6, width=1.2, labelsize=24)
        axs[i].spines["bottom"].set_linewidth(1.2)
        axs[i].spines["top"].set_linewidth(1.2)
        axs[i].spines["right"].set_linewidth(1.2)
        axs[i].spines["left"].set_linewidth(1.2)
    
# save
fig.savefig('./../fig/analysis/V_profile_{}.png'.format(file_name), format='png', dpi=300, transparent=False)
fig.savefig('./../fig/analysis/V_profile_{}.pdf'.format(file_name), format='pdf', transparent=True)

# close
plt.close()

# print
print('[Message:  5/16] Plot V_profile        has done.')

#---------------------
# V wall
#---------------------
# figure and axis environment
fig, axs = plt.subplots(2, 2, figsize=(13, 12), facecolor='white', subplot_kw={'facecolor':'white'})
# margin between figures
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.1, top=0.95, wspace=0.3, hspace=0.2)

axs[0,0].scatter(MP_b[:,0,-1]/length, np.abs(MP_b[:,2,-1])/V_nor, s=26, c='red',  label='bottom')
axs[0,0].scatter(MP_t[:,0,-1]/length, np.abs(MP_t[:,2,-1])/V_nor, s=26, c='blue', label='top')

axs[1,0].scatter(MP_b[:,0,-1]/length, np.abs(MP_b[:,3,-1])/V_nor, s=26, c='red')
axs[1,0].scatter(MP_t[:,0,-1]/length, np.abs(MP_t[:,3,-1])/V_nor, s=26, c='blue')

axs[0,1].scatter(MP_r[:,1,-1]/length, np.abs(MP_r[:,2,-1])/V_nor, s=26, c='green',  label='right')
axs[0,1].scatter(MP_l[:,1,-1]/length, np.abs(MP_l[:,2,-1])/V_nor, s=26, c='orange', label='left')

axs[1,1].scatter(MP_r[:,1,-1]/length, np.abs(MP_r[:,3,-1])/V_nor, s=26, c='green',  label='right')
axs[1,1].scatter(MP_l[:,1,-1]/length, np.abs(MP_l[:,3,-1])/V_nor, s=26, c='orange', label='left')

# axis labels
axs[1,0].set_xlabel(r'$x^{*}$',   fontsize=28, labelpad=20)
axs[1,1].set_xlabel(r'$y^{*}$',   fontsize=28, labelpad=20)
axs[0,0].set_ylabel(r'$|u^{*}|$', fontsize=28, labelpad=22)
axs[1,0].set_ylabel(r'$|v^{*}|$', fontsize=28, labelpad=22)

# lim
axs[0,1].set_xlim(0-0.05, 1+0.05)
axs[1,1].set_xlim(0-0.05, 1+0.05)

# legend
axs[0,0].legend(fontsize=18, fancybox=True, edgecolor='silver')
axs[0,1].legend(fontsize=18, fancybox=True, edgecolor='silver')


# log scale, grid, ticks, outer_frame, lines
for i in range(2):
        for j in range(2):
                axs[i,j].set_yscale('log')
                axs[i,j].grid(which='major', color='silver', linewidth=1)
                axs[i,j].tick_params(axis='both', which='major', direction='out', length=10, width=1.2, labelsize=22)
                axs[i,j].tick_params(axis='both', which='minor', direction='out', length=6,  width=1, labelsize=22)
                axs[i,j].spines["bottom"].set_linewidth(1.2)
                axs[i,j].spines["top"].set_linewidth(1.2)
                axs[i,j].spines["right"].set_linewidth(1.2)
                axs[i,j].spines["left"].set_linewidth(1.2)
                #axs[i,j].axvline(x=0, ymin=0, ymax=1, color='silver', linewidth=0.8, linestyle='--')

# title
axs[0,0].set_title(r'bottom \& top', fontsize=26, color='k', y=1.03)
axs[0,1].set_title(r'right \& left', fontsize=26, color='k', y=1.03)

# save
fig.savefig('./../fig/analysis/V_wall_{}.png'.format(file_name), format='png', dpi=300, transparent=False)
fig.savefig('./../fig/analysis/V_wall_{}.pdf'.format(file_name), format='pdf', transparent=True)

# close
plt.close()

# print
print('[Message:  6/16] Plot V_wall           has done.')

#---------------------
# local Nu fin
#---------------------
Nu_bot_ave = np.sum(np.abs(MP_b[:,6,-1]))/len(MP_b[:,6,-1])
Nu_top_ave = np.sum(np.abs(MP_t[:,6,-1]))/len(MP_t[:,6,-1])

# figure and axis environment
fig, axs = plt.subplots(1, 2, figsize=(16, 7.3), facecolor='white', subplot_kw={'facecolor':'white'})

# margin between figures
plt.subplots_adjust(left=0.12, right=0.95, bottom=0.15, top=0.9, wspace=0.5, hspace=0.3)

# plotting
axs[0].plot(MP_b[:,0,-1]/length, np.abs(MP_b[:,6,-1]), linestyle='-', marker='o', c='k')
axs[1].plot(MP_t[:,0,-1]/length, np.abs(MP_t[:,6,-1]), linestyle='-', marker='o', c='k')

# axis labels
axs[0].set_xlabel(r'$x^{*}$',   fontsize=30, labelpad=20)
axs[1].set_xlabel(r'$x^{*}$',   fontsize=30, labelpad=20)
axs[0].set_ylabel(r'$Nu$',      fontsize=30, labelpad=20)

# title
axs[0].set_title(r'$Nu_{}=~{:.5f}$'.format('h', Nu_bot_ave), fontsize=24, color='k', y=1.03)
axs[1].set_title(r'$Nu_{}=~{:.5f}$'.format('c', Nu_top_ave), fontsize=24, color='k', y=1.03)

# grid, ticks, outer_frame, lines
for i in range(2):
        axs[i].grid(which='major', color='silver', linewidth=1)
        axs[i].tick_params(axis='both', which='major', direction='out', length=10, width=1.2, labelsize=22)
        axs[i].spines["bottom"].set_linewidth(1.2)
        axs[i].spines["top"].set_linewidth(1.2)
        axs[i].spines["right"].set_linewidth(1.2)
        axs[i].spines["left"].set_linewidth(1.2)
        axs[i].set_xlim([0, 1])

# save
fig.savefig('./../fig/analysis/local_Nu_{}.png'.format(file_name), format='png', dpi=300, transparent=False)
fig.savefig('./../fig/analysis/local_Nu_{}.pdf'.format(file_name), format='pdf', transparent=True)

# close
plt.close()

# print
print('[Message:  7/16] Plot local_Nu         has done.')

#---------------------
# Nu_ave time-change
#---------------------
Nu_bot_time = np.zeros((num_save_file+1)) # not including initial value
Nu_top_time = np.zeros((num_save_file+1))
for i in range(1,num_save_file+1):
        Nu_bot_time[i] = np.sum(np.abs(MP_b[:,6,i]))/len(MP_b[:,6,i])
        Nu_top_time[i] = np.sum(np.abs(MP_t[:,6,i]))/len(MP_t[:,6,i])
# not including initial value
Nu_bot_time = np.delete(Nu_bot_time, 0, axis=0)
Nu_top_time = np.delete(Nu_top_time, 0, axis=0)

# figure and axis environment
fig, axs = plt.subplots(1, 2, figsize=(16, 7.3), facecolor='white', subplot_kw={'facecolor':'white'})

# margin between figures
plt.subplots_adjust(left=0.12, right=0.95, bottom=0.15, top=0.9, wspace=0.5, hspace=0.3)

# plotting
axs[0].plot(pass_time[1:], Nu_bot_time[:], linestyle='-', linewidth=2.5, c='k')
axs[1].plot(pass_time[1:], Nu_top_time[:], linestyle='-', linewidth=2.5, c='k')

# axis labels
axs[0].set_xlabel(r'time $(\mathrm{s})$',    fontsize=26, labelpad=18)
axs[1].set_xlabel(r'time $(\mathrm{s})$',    fontsize=26, labelpad=18)
axs[0].set_ylabel(r'$Nu$',                   fontsize=26, labelpad=20)

# title
axs[0].set_title(r'$Nu_{}=~{:.5f}$'.format('h', Nu_bot_ave), fontsize=24, color='k', y=1.03)
axs[1].set_title(r'$Nu_{}=~{:.5f}$'.format('c', Nu_top_ave), fontsize=24, color='k', y=1.03)

# grid, ticks, outer_frame, lines
for i in range(2):
        axs[i].grid(which='major', color='silver', linewidth=1)
        axs[i].tick_params(axis='both', which='major', direction='out', length=10, width=1.2, labelsize=22)
        axs[i].tick_params(axis='both', which='minor', direction='out', length=10, width=1.2, labelsize=22)
        axs[i].spines["bottom"].set_linewidth(1.2)
        axs[i].spines["top"].set_linewidth(1.2)
        axs[i].spines["right"].set_linewidth(1.2)
        axs[i].spines["left"].set_linewidth(1.2)

# save
fig.savefig('./../fig/analysis/Nu_time_change_{}.png'.format(file_name), format='png', dpi=300, transparent=False)
fig.savefig('./../fig/analysis/Nu_time_change_{}.pdf'.format(file_name), format='pdf', transparent=True)

# close
plt.close()

# print
print('[Message:  8/16] Plot Nu_time_change   has done.')

#---------------------
# error rho
#---------------------
error_rho_ave = np.zeros((num_save_file+1))
error_par_con = np.zeros((num_save_file+1))
for i in range(num_save_file+1):
        for j in range(total_num):
                tmp_rho_ave = (np.abs(SP_inner[j, 4, i] - SP_inner[j, 8, i]) / rho_0) / total_num
                tmp_par_con = (np.abs(SP_inner[j, 10, i] - 1.0) / 1.0) / total_num
                error_rho_ave[i] += tmp_rho_ave
                error_par_con[i] += tmp_par_con

# figure and axis environment
fig, axs = plt.subplots(1, 2, figsize=(16, 7.3), facecolor='white', subplot_kw={'facecolor':'white'})

# margin between figures
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.9, wspace=0.5, hspace=0.3)

# plotting
axs[0].plot(pass_time, error_rho_ave*100, color='k', linewidth=2.5, linestyle='-' )
axs[1].plot(pass_time, error_par_con*100, color='k', linewidth=2.5, linestyle='-' )

# axis labels
axs[0].set_xlabel(r'time $(\mathrm{s})$',                fontsize=26, labelpad=18)
axs[0].set_ylabel(r'$\varepsilon(\rho^{S}, \bar{\rho})~\%$', fontsize=26, labelpad=18)
axs[1].set_xlabel(r'time $(\mathrm{s})$',                fontsize=26, labelpad=18)
axs[1].set_ylabel(r'$\varepsilon(\varpi, 1)~\%$',              fontsize=26, labelpad=18)

axs[0].set_title(r'$\varepsilon_{}={:.3f}~(\%)$'.format('{\mathrm{fin}}', error_rho_ave[-1]*100), fontsize=26, color='k', y=1.03)
axs[1].set_title(r'$\varepsilon_{}={:.3f}~(\%)$'.format('{\mathrm{fin}}', error_par_con[-1]*100), fontsize=26, color='k', y=1.03)

# grid, ticks, outer_frame
for i in range(2):
        axs[i].grid(which='major', color='silver', linewidth=1)
        axs[i].tick_params(axis='both', which='major', direction='out', length=10, width=1.2, labelsize=24)
        axs[i].tick_params(axis='both', which='minor', direction='out', length=10, width=1.2, labelsize=24)
        axs[i].spines["bottom"].set_linewidth(1.2)
        axs[i].spines["top"].set_linewidth(1.2)
        axs[i].spines["right"].set_linewidth(1.2)
        axs[i].spines["left"].set_linewidth(1.2)

# save
fig.savefig('./../fig/analysis/error_rho_{}.png'.format(file_name), format='png', dpi=300, transparent=False)
fig.savefig('./../fig/analysis/error_rho_{}.pdf'.format(file_name), format='pdf', transparent=True)

# close
plt.close()

# print
print('[Message:  9/16] Plot error_rho        has done.')

#---------------------
# tem profile
#---------------------
# figure and axis environment
fig, axs = plt.subplots(1, 2, figsize=(16, 7.3), facecolor='white', subplot_kw={'facecolor':'white'})

# margin between figures
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.18, top=0.93, wspace=0.5, hspace=0.3)

# plotting
axs[0].plot((MP_vm[:,5,-1]-T_0)/delta_tem, MP_vm[:,1,-1]/length, color='k', linewidth=2.5, linestyle='-' )
axs[1].plot((MP_r[:,5,-1] -T_0)/delta_tem, MP_r[:,1,-1]/length,  color='r', linewidth=2.5, linestyle='-' , label='right')
axs[1].plot((MP_l[:,5,-1] -T_0)/delta_tem, MP_l[:,1,-1]/length,  color='b', linewidth=2.5, linestyle='-' , label='left')

axs[1].legend(fontsize=18, fancybox=True, edgecolor='silver')
# grid, ticks, outer_frame
for i in range(2):
        axs[i].grid(which='major', color='silver', linewidth=1)
        axs[i].tick_params(axis='both', which='major', direction='out', length=10, width=1.2, labelsize=24)
        axs[i].spines["bottom"].set_linewidth(1.2)
        axs[i].spines["top"].set_linewidth(1.2)
        axs[i].spines["right"].set_linewidth(1.2)
        axs[i].spines["left"].set_linewidth(1.2)
        axs[i].set_xlabel(r'$T^{*}$', fontsize=30, labelpad=18)
        axs[i].set_ylabel(r'$y^{*}$', fontsize=32, labelpad=22)

# save
fig.savefig('./../fig/analysis/tem_profile_{}.png'.format(file_name), format='png', dpi=300, transparent=False)
fig.savefig('./../fig/analysis/tem_profile_{}.pdf'.format(file_name), format='pdf', transparent=True)

# close
plt.close()

# print
print('[Message: 10/16] Plot tem_profile      has done.')

#---------------------
# final PS_shift
#---------------------
i = num_save_file
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
# axis labels, ticks, outer frame, lines
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
# title
plt.suptitle('{:12.3e}s  ({}step)'.format(pass_time[i], eff_step), fontsize=30, color='k', x=0.515, y=0.975)

# save
fig.savefig('./../fig/analysis/PS_shift_{}.png'.format(file_name), format='png', dpi=300, transparent=False)
fig.savefig('./../fig/analysis/PS_shift_{}.pdf'.format(file_name), format='pdf', transparent=True)

# close
plt.close()
print('[Message: 11/16] Plot_PS_shift         has done.')

# #---------------------
# # RSS change
# #---------------------
# # figure and axis environment
# fig, axs = plt.subplots(1, 2, figsize=(16, 7.3), facecolor='white', subplot_kw={'facecolor':'white'})

# # margin between figures
# plt.subplots_adjust(left=0.15, right=0.95, bottom=0.18, top=0.92, wspace=0.5, hspace=0.3)

# tmp_arg = np.zeros((time_step.shape[0]))
# for i in range(1, len(tmp_arg)):
#         tmp_arg[i] += i * write_step

# # plotting
# axs[0].plot(tmp_arg, time_step, color='k', linewidth=2.5, linestyle='-' )
# axs[1].plot(tmp_arg, RSS_arg,   color='k', linewidth=2.5, linestyle='-' , label='right')

# # label
# axs[0].set_xlabel('step', fontsize=30, labelpad=18)
# axs[0].set_ylabel(r'$\Delta t$', fontsize=32, labelpad=22)
# axs[1].set_xlabel('step', fontsize=30, labelpad=18)
# axs[1].set_ylabel(r'$c_{\zeta\xi}$', fontsize=32, labelpad=22)

# # grid, ticks, outer_frame
# for i in range(2):
#         axs[i].set_yscale('log')
#         axs[i].grid(which='major', color='silver', linewidth=1)
#         axs[i].tick_params(axis='both', which='major', direction='out', length=10, width=1.2, labelsize=24)
#         axs[i].tick_params(axis='both', which='minor', direction='out', length=10, width=1.2, labelsize=24)
#         axs[i].spines["bottom"].set_linewidth(1.2)
#         axs[i].spines["top"].set_linewidth(1.2)
#         axs[i].spines["right"].set_linewidth(1.2)
#         axs[i].spines["left"].set_linewidth(1.2)

# # save
# fig.savefig('./../fig/analysis/RSS_change_{}.png'.format(file_name), format='png', dpi=300, transparent=False)
# fig.savefig('./../fig/analysis/RSS_change_{}.pdf'.format(file_name), format='pdf', transparent=True)

# # close
# plt.close()

# # print
# print('[Message: 12/16] Plot RSS_change       has done.')

#---------------------
# Re scatter
#---------------------
i = num_save_file
eff_step = i * write_step

# figure and axis environment
fig, axs = plt.subplots(1, 2, figsize=(15, 9), facecolor='white', subplot_kw={'facecolor':'white'})
        
# margin between figure
plt.subplots_adjust(left=0.1, right=0.95, bottom=0.02, top=0.95, wspace=0.35, hspace=0.1)
        
# Re scatter
fig_1 = axs[0].scatter(SP_inner[:,0,i]/length, SP_inner[:,1,i]/length, s=SP_size, \
        c=np.abs(EOM_inner[:,5,i]/EOM_inner[:,0,i]), cmap='jet')
bar_1 = plt.colorbar(fig_1, aspect=60, ax=axs[0], orientation='horizontal', pad=0.18)
bar_1.set_label(r'$Re_{\mathrm{eff}}~(x)$', size=28, labelpad=14)
bar_1.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=24)

fig_2 = axs[1].scatter(SP_inner[:,0,i]/length, SP_inner[:,1,i]/length, s=SP_size, \
        c=np.abs(EOM_inner[:,6,i]/EOM_inner[:,1,i]), cmap='jet')
bar_2 = plt.colorbar(fig_2, aspect=60, ax=axs[1], orientation='horizontal', pad=0.18)
bar_2.set_label(r'$Re_{\mathrm{eff}}~(y)$', size=28, labelpad=14)
bar_2.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=24)

for i in range(2):
        axs[i].set_xlabel(r'$x^{*}$', fontsize=30, labelpad=14)
        axs[i].set_ylabel(r'$y^{*}$', fontsize=30, labelpad=20)

        axs[i].axvline(x=0, ymin=0, ymax=1, color='silver', linewidth=0.8, linestyle='--')
        axs[i].axvline(x=1, ymin=0, ymax=1, color='silver', linewidth=0.8, linestyle='--')
        axs[i].axhline(y=0, xmin=0, xmax=1, color='silver', linewidth=0.8, linestyle='--')
        axs[i].axhline(y=1, xmin=0, xmax=1, color='silver', linewidth=0.8, linestyle='--')

        axs[i].tick_params(axis='both', which='major', direction='out', length=3, width=0.8, labelsize=24)
        axs[i].tick_params(axis='both', which='minor', direction='out', length=3, width=0.8, labelsize=24)
        axs[i].spines["bottom"].set_linewidth(1.2)
        axs[i].spines["top"].set_linewidth(1.2)
        axs[i].spines["right"].set_linewidth(1.2)
        axs[i].spines["left"].set_linewidth(1.2)

# save
fig.savefig('./../fig/analysis/Re_scatter_{}.png'.format(file_name), format='png', dpi=300, transparent=False)
fig.savefig('./../fig/analysis/Re_scatter_{}.pdf'.format(file_name), format='pdf', transparent=True)

# close
plt.close()
print('[Message: 13/16] Plot_Re_scatter       has done.')

#---------------------
# right hand of EOM
#---------------------
# not including initial state
vis_EOM = np.zeros((num_save_file, 2))
buo_EOM = np.zeros((num_save_file))
pre_EOM = np.zeros((num_save_file, 2))
ine_EOM = np.zeros((num_save_file, 2))

for i in range(num_save_file):
        vis_EOM[i,0] = np.sum(np.abs(EOM_inner[:,0,i+1]))
        vis_EOM[i,1] = np.sum(np.abs(EOM_inner[:,1,i+1]))
        buo_EOM[i]   = np.sum(np.abs(EOM_inner[:,2,i+1]))
        pre_EOM[i,0] = np.sum(np.abs(EOM_inner[:,3,i+1]))
        pre_EOM[i,1] = np.sum(np.abs(EOM_inner[:,4,i+1]))
        ine_EOM[i,0] = np.sum(np.abs(EOM_inner[:,5,i+1]))
        ine_EOM[i,1] = np.sum(np.abs(EOM_inner[:,6,i+1]))

# figure and axis environment
fig, axs = plt.subplots(1, 2, figsize=(16, 8), facecolor='white', subplot_kw={'facecolor':'white'})
        
# margin between figure
plt.subplots_adjust(left=0.1, right=0.95, bottom=0.15, top=0.9, wspace=0.4, hspace=0.1)

axs[0].plot(pass_time[1:], vis_EOM[:,0], color='k', linewidth=2.5, linestyle='-', label='vis')
axs[0].plot(pass_time[1:], pre_EOM[:,0], color='b', linewidth=2.5, linestyle='-', label='pre')
axs[0].plot(pass_time[1:], ine_EOM[:,0], color='r', linewidth=2.5, linestyle='-', label='ine')

axs[1].plot(pass_time[1:], vis_EOM[:,1], color='k', linewidth=2.5, linestyle='-', label='vis')
axs[1].plot(pass_time[1:], pre_EOM[:,1], color='b', linewidth=2.5, linestyle='-', label='pre')
axs[1].plot(pass_time[1:], ine_EOM[:,1], color='r', linewidth=2.5, linestyle='-', label='ine')
axs[1].plot(pass_time[1:], buo_EOM[:],   color='magenta', linewidth=2.5, linestyle='-', label='buo')

for i in range(2):
        axs[i].grid(which='major', color='silver', linewidth=1)
        axs[i].tick_params(axis='both', which='major', direction='out', length=3, width=0.8, labelsize=24)
        axs[i].tick_params(axis='both', which='minor', direction='out', length=3, width=0.8, labelsize=24)
        axs[i].spines["bottom"].set_linewidth(1.2)
        axs[i].spines["top"].set_linewidth(1.2)
        axs[i].spines["right"].set_linewidth(1.2)
        axs[i].spines["left"].set_linewidth(1.2)
        axs[i].set_xlabel('time (s)', fontsize=24, labelpad=14)
        axs[i].set_xlabel(r'time $(\mathrm{s})$', fontsize=26, labelpad=18)
        axs[i].legend(fontsize=24, fancybox=True, edgecolor='silver')

axs[0].set_ylabel(r'$\mathrm{EOM}_{x}~(\mathrm{m~s^{-2}})$', fontsize=26, labelpad=20)
axs[1].set_ylabel(r'$\mathrm{EOM}_{y}~(\mathrm{m~s^{-2}})$', fontsize=26, labelpad=20)

# save
fig.savefig('./../fig/analysis/EOM_{}.png'.format(file_name), format='png', dpi=300, transparent=False)
fig.savefig('./../fig/analysis/EOM_{}.pdf'.format(file_name), format='pdf', transparent=True)

# close
plt.close()
print('[Message: 14/16] Plot_EOM              has done.')

#---------------------
# Re_change
#---------------------
Re_tmp = np.zeros((total_num, 2, num_save_file))
Re_max = np.zeros((num_save_file, 2))

for i in range(num_save_file):
        Re_tmp[:,0,i] = np.abs(EOM_inner[:,5,i+1]/EOM_inner[:,0,i+1])
        Re_tmp[:,1,i] = np.abs(EOM_inner[:,6,i+1]/EOM_inner[:,1,i+1])

for i in range(num_save_file):
        Re_max[i,0] = np.max(Re_tmp[:,0,i], axis=0)
        Re_max[i,1] = np.max(Re_tmp[:,1,i], axis=0)

# figure and axis environment
fig, axs = plt.subplots(1, 2, figsize=(15, 7.5), facecolor='white', subplot_kw={'facecolor':'white'})
        
# margin between figure
plt.subplots_adjust(left=0.10, right=0.95, bottom=0.15, top=0.9, wspace=0.35, hspace=0.1)

axs[0].plot(pass_time[1:], Re_max[:,0], color='k', linewidth=2.5, linestyle='-')
axs[1].plot(pass_time[1:], Re_max[:,1], color='k', linewidth=2.5, linestyle='-')

for i in range(2):
        axs[i].grid(which='major', color='silver', linewidth=1)
        axs[i].tick_params(axis='both', which='major', direction='out', length=3, width=0.8, labelsize=24)
        axs[i].tick_params(axis='both', which='minor', direction='out', length=3, width=0.8, labelsize=24)
        axs[i].spines["bottom"].set_linewidth(1.2)
        axs[i].spines["top"].set_linewidth(1.2)
        axs[i].spines["right"].set_linewidth(1.2)
        axs[i].spines["left"].set_linewidth(1.2)
        axs[i].set_xlabel('time (s)', fontsize=24, labelpad=14)
        axs[i].set_xlabel(r'time $(\mathrm{s})$', fontsize=26, labelpad=18)

axs[0].set_ylabel(r'$Re_{\mathrm{eff}}~(x)$', fontsize=26, labelpad=18)
axs[1].set_ylabel(r'$Re_{\mathrm{eff}}~(y)$', fontsize=26, labelpad=18)

# title
axs[0].set_title(r'$Re_{}={:.2e}$'.format('{\mathrm{eff, fin}}', Re_max[-1,0]), fontsize=24, color='k', y=1.02)
axs[1].set_title(r'$Re_{}={:.2e}$'.format('{\mathrm{eff, fin}}', Re_max[-1,1]), fontsize=24, color='k', y=1.02)

# save
fig.savefig('./../fig/analysis/Re_change_{}.png'.format(file_name), format='png', dpi=300, transparent=False)
fig.savefig('./../fig/analysis/Re_change_{}.pdf'.format(file_name), format='pdf', transparent=True)

# close
plt.close()
print('[Message: 15/16] Plot_Re_change        has done.')

print('[Message: 16/16] plot_analysis         has completed !')
end_time = time.perf_counter()
print('-----------------------------------------------------')
print('[Message] Total time : {:.2f} [s]'.format(end_time-start_time))
print('-----------------------------------------------------')
print('')