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
print('[Message]  plot_initial_state has started.')
print('-----------------------------------------------------')
start_time = time.perf_counter()

#=====================
# scatter size
#=====================
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

#=====================
# read data
#=====================
# get file_name
file_path = glob.glob('./../output_setting/initial_state_SP_*') # finding file
file_name = file_path[0][-8:-4] # [-8:-4] is file_name in str type
endian = '>'

# SP
f = open('./../output_setting/initial_state_SP_{}.dat'.format(file_name), 'rb')
data_type = np.dtype([('SP_xy',     endian+str(2*total_num_system)+'d'), \
                      ('SP_uv',     endian+str(2*total_num_system)+'d'), \
                      ('SP_rho_S',  endian+str(total_num_system)+'d'), \
                      ('SP_pre',    endian+str(total_num_system)+'d'), \
                      ('SP_tem',    endian+str(total_num_system)+'d'), \
                      ('SP_eta',    endian+str(total_num_system)+'d'), \
                      ('SP_sm_rho', endian+str(total_num_system)+'d'), \
                      ('SP_rho_G',  endian+str(total_num_system)+'d')])
data = np.fromfile(f, dtype=data_type)
SP_xy     = data['SP_xy'].reshape(    (total_num_system, 2), order='F')
SP_uv     = data['SP_uv'].reshape(    (total_num_system, 2), order='F')
SP_rho_S  = data['SP_rho_S'].reshape( (total_num_system, 1), order='F')
SP_pre    = data['SP_pre'].reshape(   (total_num_system, 1), order='F')
SP_tem    = data['SP_tem'].reshape(   (total_num_system, 1), order='F')
SP_eta    = data['SP_eta'].reshape(   (total_num_system, 1), order='F')
SP_sm_rho = data['SP_sm_rho'].reshape((total_num_system, 1), order='F')
SP_rho_G  = data['SP_rho_G'].reshape( (total_num_system, 1), order='F')

# particle ID
f = open('./../output_setting/initial_state_SP_kind_{}.dat'.format(file_name), 'rb')
data_type = np.dtype([('SP_kind', endian+str(total_num_system)+'i')])
data = np.fromfile(f, dtype=data_type)
SP_kind = data['SP_kind'].reshape((total_num_system, 1), order='F')

# NNPS
f = open('./../output_setting/initial_state_NNPS_{}.dat'.format(file_name), 'rb')
data_type = np.dtype([('SP_uv', endian+str(total_num_system)+'i')])
data = np.fromfile(f, dtype=data_type)
NNPS = data['SP_uv'].reshape((total_num_system, 1), order='F')

# VM
f = open('./../output_setting/initial_VM_{}.dat'.format(file_name), 'rb')
data_type = np.dtype([('VM_xy',    endian+str(2*num_wall)+'d'), \
                      ('VM_uv',    endian+str(2*num_wall)+'d'), \
                      ('VM_pre',   endian+str(num_wall)+'d'), \
                      ('VM_tem',   endian+str(num_wall)+'d'), \
                      ('ID_VM_WL', endian+str(num_wall)+'i')])
data = np.fromfile(f, dtype=data_type)
VM_xy    = data['VM_xy'].reshape(   (num_wall, 2), order='F')
VM_uv    = data['VM_uv'].reshape(   (num_wall, 2), order='F')
VM_pre   = data['VM_pre'].reshape(  (num_wall, 1), order='F')
VM_tem   = data['VM_tem'].reshape(  (num_wall, 1), order='F')
ID_VM_WL = data['ID_VM_WL'].reshape((num_wall, 1), order='F')

# NNPS for VM
f = open('./../output_setting/NNPS_VM_{}.dat'.format(file_name), 'rb')
data_type = np.dtype([('VM_uv', endian+str(num_wall)+'i')])
data = np.fromfile(f, dtype=data_type)
NNPS_VM = data['VM_uv'].reshape((num_wall, 1), order='F')

# MP
f = open('./../output_setting/initial_MP_{}.dat'.format(file_name), 'rb')
data_type = np.dtype([('MP_xy',  endian+str(2*num_MP)+'d'), \
                      ('MP_uv',  endian+str(2*num_MP)+'d'), \
                      ('MP_pre', endian+str(num_MP)+'d'), \
                      ('MP_tem', endian+str(num_MP)+'d'), \
                      ('MP_NU',  endian+str(num_MP)+'d') ])
data = np.fromfile(f, dtype=data_type)
MP_xy  = data['MP_xy'].reshape( (num_MP, 2), order='F')
MP_uv  = data['MP_uv'].reshape( (num_MP, 2), order='F')
MP_pre = data['MP_pre'].reshape((num_MP, 1), order='F')
MP_tem = data['MP_tem'].reshape((num_MP, 1), order='F')
MP_NU  = data['MP_NU'].reshape( (num_MP, 1), order='F')

# NNPS for MP
f = open('./../output_setting/NNPS_MP_{}.dat'.format(file_name), 'rb')
data_type = np.dtype([('MP_uv', endian+str(num_MP)+'i')])
data = np.fromfile(f, dtype=data_type)
NNPS_MP = data['MP_uv'].reshape((num_MP, 1), order='F')

# translate axes
tra_dis = wall_thickness * delta_x
SP_xy[:, 0] = SP_xy[:, 0] - tra_dis
SP_xy[:, 1] = SP_xy[:, 1] - tra_dis
VM_xy[:, 0] = VM_xy[:, 0] - tra_dis
VM_xy[:, 1] = VM_xy[:, 1] - tra_dis
MP_xy[:, 0] = MP_xy[:, 0] - tra_dis
MP_xy[:, 1] = MP_xy[:, 1] - tra_dis

# concatenate
SP_read = np.hstack((SP_xy, SP_uv, SP_rho_S, SP_pre, SP_tem, SP_eta, SP_sm_rho, SP_rho_G))
VM_read = np.hstack((VM_xy, VM_uv, VM_pre, VM_tem, ID_VM_WL))
MP_read = np.hstack((MP_xy, MP_uv, MP_pre, MP_tem))

# categorize
SP_inner = np.zeros((1, 10))
SP_wall  = np.zeros((1, 10))
VM_corner = np.zeros((1, 7))
NNPS_VM_corner = np.zeros((1))
MP_b  = np.zeros((1, 6))
MP_hm = np.zeros((1, 6))
MP_t  = np.zeros((1, 6))
MP_l  = np.zeros((1, 6))
MP_vm = np.zeros((1, 6))
MP_r  = np.zeros((1, 6))

for i in range(total_num_system):
    if (SP_kind[i] == 0):
        SP_inner = np.append(SP_inner, [SP_read[i, :]], axis=0)
    else:
        SP_wall = np.append(SP_wall, [SP_read[i, :]], axis=0)

tmp = wall_thickness * delta_x
for i in range(VM_read.shape[0]):
    you = int(VM_read[i, 6]) - 1
    if (SP_read[you, 0]/length < 0) and (SP_read[you, 1]/length < 0):
        VM_corner      = np.append(VM_corner, [VM_read[i, :]], axis=0)
        NNPS_VM_corner = np.append(NNPS_VM_corner, NNPS_VM[i], axis=0)
    elif (SP_read[you, 0]/length > 1) and (SP_read[you, 1]/length < 0):
        VM_corner      = np.append(VM_corner, [VM_read[i, :]], axis=0)
        NNPS_VM_corner = np.append(NNPS_VM_corner, NNPS_VM[i], axis=0)
    elif (SP_read[you, 0]/length < 0) and (SP_read[you, 1]/length > 1):
        VM_corner      = np.append(VM_corner, [VM_read[i, :]], axis=0)
        NNPS_VM_corner = np.append(NNPS_VM_corner, NNPS_VM[i], axis=0)
    elif (SP_read[you, 0]/length > 1) and (SP_read[you, 1]/length > 1):
        VM_corner      = np.append(VM_corner, [VM_read[i, :]], axis=0)
        NNPS_VM_corner = np.append(NNPS_VM_corner, NNPS_VM[i], axis=0)

for i in range(num_MP):
    if (0*num_MP_side <= i < 1*num_MP_side):
        MP_b = np.append(MP_b, [MP_read[i, :]], axis=0)
    elif (1*num_MP_side <= i < 2*num_MP_side):
        MP_hm = np.append(MP_hm, [MP_read[i, :]], axis=0)
    elif (2*num_MP_side <= i < 3*num_MP_side):
        MP_t = np.append(MP_t, [MP_read[i, :]], axis=0)
    elif (3*num_MP_side <= i < 4*num_MP_side):
        MP_l = np.append(MP_l, [MP_read[i, :]], axis=0)
    elif (4*num_MP_side <= i < 5*num_MP_side):
        MP_vm = np.append(MP_vm, [MP_read[i, :]], axis=0)
    elif (5*num_MP_side <= i < 6*num_MP_side):
        MP_r = np.append(MP_r, [MP_read[i, :]], axis=0)

SP_inner  = np.delete(SP_inner, 0, axis=0)
SP_wall   = np.delete(SP_wall, 0, axis=0)
VM_corner = np.delete(VM_corner, 0, axis=0)
NNPS_VM_corner = np.delete(NNPS_VM_corner, 0, axis=0)
MP_b  = np.delete(MP_b, 0, axis=0)
MP_hm = np.delete(MP_hm, 0, axis=0)
MP_t  = np.delete(MP_t, 0, axis=0)
MP_l  = np.delete(MP_l, 0, axis=0)
MP_vm = np.delete(MP_vm, 0, axis=0)
MP_r =  np.delete(MP_r, 0, axis=0)

print('[Message: 1/7] Reading data         has done.')

#=====================
# plot initial state
#=====================
# figure and axis environment
fig, axs = plt.subplots(2, 2, figsize=(12, 14), facecolor='white', subplot_kw={'facecolor':'white'})
# margin between figures
plt.subplots_adjust(left=0.12, right=0.9, bottom=0.05, top=0.95, wspace=0.4, hspace=0.1)

# plot
fig_1 = axs[0,0].scatter(SP_read[:,0]/length, SP_read[:,1]/length, s=SP_size, c=SP_kind[:], cmap='tab20b', vmin=0, vmax=4)
bar_1 = plt.colorbar(fig_1, aspect=60, ax=axs[0,0], orientation='horizontal', pad=0.16)
bar_1.set_label('particle ID \n(0:inner, 1:bottom, 2:top, 3:right, 4:left)', size=18, labelpad=12)
bar_1.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)

fig_2 = axs[1,0].scatter(SP_read[:,0]/length, SP_read[:,1]/length, s=SP_size, c=(SP_read[:,6]-T_0)/delta_tem, cmap='jet', vmin=-0.501, vmax=0.501)
bar_2 = plt.colorbar(fig_2, aspect=60, ax=axs[1,0], orientation='horizontal', pad=0.16)
bar_2.set_label(r'$T^{*}$', size=18, labelpad=12)
bar_2.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)

fig_3 = axs[0,1].scatter(SP_inner[:,0]/length, SP_inner[:,1]/length, s=SP_size, c=(SP_inner[:,9] + rho_0)/rho_0, cmap='jet')
bar_3 = plt.colorbar(fig_3, aspect=60, ax=axs[0,1], orientation='horizontal', pad=0.16)
bar_3.set_label(r'$\rho^{G*}$', size=18, labelpad=12)
bar_3.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)

fig_4 = axs[1,1].scatter(SP_read[:,0]/length, SP_read[:,1]/length, s=SP_size, c=NNPS[:], cmap='prism', vmin=0)
bar_4 = plt.colorbar(fig_4, aspect=60, ax=axs[1,1], orientation='horizontal', pad=0.16)
bar_4.set_label('NNPS number', size=18, labelpad=12)
bar_4.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)

# axs_label, tick_prams, outer_frame, lines
for i in range(2):
    for j in range(2):
        axs[i,j].set_xlabel(r'$x^{*}$', fontsize=22, labelpad=14)
        axs[i,j].set_ylabel(r'$y^{*}$', fontsize=22, labelpad=22)

        axs[i,j].tick_params(axis='both', which='major', direction='out', length=3, width=0.8, labelsize=16)

        axs[i,j].spines["bottom"].set_linewidth(1.2)
        axs[i,j].spines["top"].set_linewidth(1.2)
        axs[i,j].spines["right"].set_linewidth(1.2)
        axs[i,j].spines["left"].set_linewidth(1.2)

        axs[i,j].axvline(x=0, ymin=0, ymax=1, color='silver', linewidth=0.8, linestyle='--')
        axs[i,j].axvline(x=1, ymin=0, ymax=1, color='silver', linewidth=0.8, linestyle='--')
        axs[i,j].axhline(y=0, xmin=0, xmax=1, color='silver', linewidth=0.8, linestyle='--')
        axs[i,j].axhline(y=1, xmin=0, xmax=1, color='silver', linewidth=0.8, linestyle='--')

# save
fig.savefig('./../fig/initial_state/initial_setting_{}.png'.format(file_name), \
    format='png', dpi=300, transparent=False)
fig.savefig('./../fig/initial_state/initial_setting_{}.pdf'.format(file_name), \
    format='pdf', transparent=True)
# close
plt.close()
print('[Message: 2/7] plot_initial_state   has done.')

#=====================
# plot initial VM
#=====================
# figure and axis environment
fig, axs = plt.subplots(2, 2, figsize=(12, 14), facecolor='white', subplot_kw={'facecolor':'white'})
# margin between figures
plt.subplots_adjust(left=0.12, right=0.9, bottom=0.05, top=0.95, wspace=0.4, hspace=0.1)

# plot
tmp = wall_thickness * num_system
VM_count = VM_read.shape[0]
num_VM_sides = VM_count - tmp
fig_1 = axs[0,0].scatter(VM_read[:tmp,0]/length, VM_read[:tmp,1]/length, s=SP_size, c=VM_read[:tmp,6], cmap='tab20b', vmin=0)
bar_1 = plt.colorbar(fig_1, aspect=60, ax=axs[0,0], orientation='horizontal', pad=0.16)
bar_1.set_label('particle number', size=18, labelpad=12)
bar_1.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)

fig_2 = axs[1,0].scatter(VM_read[tmp:num_VM_sides,0]/length, VM_read[tmp:num_VM_sides,1]/length, \
    s=SP_size, c=VM_read[tmp:num_VM_sides,6], cmap='tab20b')
bar_2 = plt.colorbar(fig_2, aspect=60, ax=axs[1,0], orientation='horizontal', pad=0.16)
bar_2.set_label('particle number', size=18, labelpad=12)
bar_2.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)

fig_3 = axs[0,1].scatter(VM_read[num_VM_sides+1:,0]/length, VM_read[num_VM_sides+1:,1]/length, \
    s=SP_size, c=VM_read[num_VM_sides+1:,6], cmap='tab20b')
bar_3 = plt.colorbar(fig_3, aspect=60, ax=axs[0,1], orientation='horizontal', pad=0.16)
bar_3.set_label('particle number', size=18, labelpad=12)
bar_3.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)

fig_4 = axs[1,1].scatter(VM_corner[:,0]/length, VM_corner[:,1]/length, s=SP_size, c=VM_corner[:,6], cmap='tab20b', vmin=0)
bar_4 = plt.colorbar(fig_4, aspect=60, ax=axs[1,1], orientation='horizontal', pad=0.16)
bar_4.set_label('particle number', size=18, labelpad=12)
bar_4.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)

# axs_label, tick_prams, outer_frame, lines
for i in range(2):
    for j in range(2):
        axs[i,j].set_xlabel(r'$x^{*}$', fontsize=22, labelpad=14)
        axs[i,j].set_ylabel(r'$y^{*}$', fontsize=22, labelpad=22)

        axs[i,j].tick_params(axis='both', which='major', direction='out', length=3, width=0.8, labelsize=16)

        axs[i,j].spines["bottom"].set_linewidth(1.2)
        axs[i,j].spines["top"].set_linewidth(1.2)
        axs[i,j].spines["right"].set_linewidth(1.2)
        axs[i,j].spines["left"].set_linewidth(1.2)

        axs[i,j].axvline(x=0, ymin=0, ymax=1, color='silver', linewidth=0.8, linestyle='--')
        axs[i,j].axvline(x=1, ymin=0, ymax=1, color='silver', linewidth=0.8, linestyle='--')
        axs[i,j].axhline(y=0, xmin=0, xmax=1, color='silver', linewidth=0.8, linestyle='--')
        axs[i,j].axhline(y=1, xmin=0, xmax=1, color='silver', linewidth=0.8, linestyle='--')

# save
fig.savefig('./../fig/initial_state/initial_VM_{}.png'.format(file_name), \
    format='png', dpi=300, transparent=False)
fig.savefig('./../fig/initial_state/initial_VM_{}.pdf'.format(file_name), \
    format='pdf', transparent=True)
# close
plt.close()
print('[Message: 3/7] plot_initial_VM      has done.')

#=====================
# plot NNPS_VM
#=====================
# figure and axis environment
fig, axs = plt.subplots(2, 2, figsize=(12, 14), facecolor='white', subplot_kw={'facecolor':'white'})
# margin between figures
plt.subplots_adjust(left=0.12, right=0.9, bottom=0.1, top=0.95, wspace=0.4, hspace=0.1)

# plot
tmp = wall_thickness * num_system
VM_count = VM_read.shape[0]
num_VM_sides = VM_count - tmp
fig_1 = axs[0,0].scatter(VM_read[:tmp,0]/length, VM_read[:tmp,1]/length, s=SP_size, c=NNPS_VM[:tmp], cmap='tab20b')
bar_1 = plt.colorbar(fig_1, aspect=60, ax=axs[0,0], orientation='horizontal', pad=0.16)
bar_1.set_label('NNPS number', size=18, labelpad=12)
bar_1.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)

fig_2 = axs[1,0].scatter(VM_read[tmp:num_VM_sides,0]/length, VM_read[tmp:num_VM_sides,1]/length, \
    s=SP_size, c=NNPS_VM[tmp:num_VM_sides], cmap='tab20b')
bar_2 = plt.colorbar(fig_2, aspect=60, ax=axs[1,0], orientation='horizontal', pad=0.16)
bar_2.set_label('NNPS number', size=18, labelpad=12)
bar_2.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)

fig_3 = axs[0,1].scatter(VM_read[num_VM_sides+1:,0]/length, VM_read[num_VM_sides+1:,1]/length, \
    s=SP_size, c=NNPS_VM[num_VM_sides+1:], cmap='tab20b')
bar_3 = plt.colorbar(fig_3, aspect=60, ax=axs[0,1], orientation='horizontal', pad=0.16)
bar_3.set_label('NNPS number', size=18, labelpad=12)
bar_3.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)

fig_4 = axs[1,1].scatter(VM_corner[:,0]/length, VM_corner[:,1]/length, s=SP_size, c=NNPS_VM_corner[:], cmap='tab20b')
bar_4 = plt.colorbar(fig_4, aspect=60, ax=axs[1,1], orientation='horizontal', pad=0.16)
bar_4.set_label('NNPS number', size=18, labelpad=12)
bar_4.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)

# axs_label, tick_prams, outer_frame, lines
for i in range(2):
    for j in range(2):
        axs[i,j].set_xlabel(r'$x^{*}$', fontsize=22, labelpad=14)
        axs[i,j].set_ylabel(r'$y^{*}$', fontsize=22, labelpad=22)

        axs[i,j].tick_params(axis='both', which='major', direction='out', length=3, width=0.8, labelsize=16)

        axs[i,j].spines["bottom"].set_linewidth(1.2)
        axs[i,j].spines["top"].set_linewidth(1.2)
        axs[i,j].spines["right"].set_linewidth(1.2)
        axs[i,j].spines["left"].set_linewidth(1.2)

        axs[i,j].axvline(x=0, ymin=0, ymax=1, color='silver', linewidth=0.8, linestyle='--')
        axs[i,j].axvline(x=1, ymin=0, ymax=1, color='silver', linewidth=0.8, linestyle='--')
        axs[i,j].axhline(y=0, xmin=0, xmax=1, color='silver', linewidth=0.8, linestyle='--')
        axs[i,j].axhline(y=1, xmin=0, xmax=1, color='silver', linewidth=0.8, linestyle='--')
# save
fig.savefig('./../fig/initial_state/NNPS_VM_{}.png'.format(file_name), \
    format='png', dpi=300, transparent=False)
fig.savefig('./../fig/initial_state/NNPS_VM_{}.pdf'.format(file_name), \
    format='pdf', transparent=True)
# close
plt.close()
print('[Message: 4/7] plot_NNPS_VM         has done.')

#=====================
# plot MP
#=====================
# figure and axis environment
fig, axs = plt.subplots(1, 3, figsize=(12, 4), facecolor='white', subplot_kw={'facecolor':'white'})
# margin between figures
plt.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.95, wspace=0.4, hspace=0.3)

# plot
axs[0].scatter(MP_b[:,0]/length, MP_b[:,1]/length, s=SP_size, c='red')
axs[0].scatter(MP_t[:,0]/length, MP_t[:,1]/length, s=SP_size, c='blue')

axs[1].scatter(MP_hm[:,0]/length, MP_hm[:,1]/length, s=SP_size, c='red')
axs[1].scatter(MP_vm[:,0]/length, MP_vm[:,1]/length, s=SP_size, c='blue')

axs[2].scatter(MP_l[:,0]/length, MP_l[:,1]/length, s=SP_size, c='red')
axs[2].scatter(MP_r[:,0]/length, MP_r[:,1]/length, s=SP_size, c='blue')

axs[0].set_ylabel(r'$y^{*}$', fontsize=22, labelpad=22)

# axs_label, tick_prams, outer_frame, lines
for i in range(3):
    axs[i].set_xlabel(r'$x^{*}$', fontsize=22, labelpad=14)

    axs[i].tick_params(axis='both', which='major', direction='out', length=3, width=0.8, labelsize=16)

    axs[i].spines["bottom"].set_linewidth(1.2)
    axs[i].spines["top"].set_linewidth(1.2)
    axs[i].spines["right"].set_linewidth(1.2)
    axs[i].spines["left"].set_linewidth(1.2)

    axs[i].axvline(x=0, ymin=0, ymax=1, color='silver', linewidth=0.8, linestyle='--')
    axs[i].axvline(x=1, ymin=0, ymax=1, color='silver', linewidth=0.8, linestyle='--')
    axs[i].axhline(y=0, xmin=0, xmax=1, color='silver', linewidth=0.8, linestyle='--')
    axs[i].axhline(y=1, xmin=0, xmax=1, color='silver', linewidth=0.8, linestyle='--')
# save
fig.savefig('./../fig/initial_state/MP_setting_{}.png'.format(file_name), \
    format='png', dpi=300, transparent=False)
fig.savefig('./../fig/initial_state/MP_setting_{}.pdf'.format(file_name), \
    format='pdf', transparent=True)
# close
plt.close()
print('[Message: 5/7] plot_MP_setting      has done.')

#=====================
# plot NNPS_MP
#=====================
# figure and axis environment
fig, axs = plt.subplots(1, 1, figsize=(6,7), facecolor='white', subplot_kw={'facecolor':'white'})
# margin between figures
plt.subplots_adjust(left=0.2, right=0.9, bottom=0.05, top=0.95, wspace=0.4, hspace=0.1)

# plot
fig_4 = axs.scatter(MP_read[:,0]/length, MP_read[:,1]/length, s=SP_size, c=MP_read[:,2], cmap='tab20b')
bar_4 = plt.colorbar(fig_4, aspect=60, ax=axs, orientation='horizontal', pad=0.16)
bar_4.set_label('NNPS number', size=18, labelpad=12)
bar_4.ax.tick_params(direction='out', length=2.5, width=0.8, labelsize=16)

# axs_label, tick_prams, outer_frame, lines
axs.set_xlabel(r'$x^{*}$', fontsize=22, labelpad=14)
axs.set_ylabel(r'$y^{*}$', fontsize=22, labelpad=22)

axs.tick_params(axis='both', which='major', direction='out', length=3, width=0.8, labelsize=16)

axs.spines["bottom"].set_linewidth(1.2)
axs.spines["top"].set_linewidth(1.2)
axs.spines["right"].set_linewidth(1.2)
axs.spines["left"].set_linewidth(1.2)

axs.axvline(x=0, ymin=0, ymax=1, color='silver', linewidth=0.8, linestyle='--')
axs.axvline(x=1, ymin=0, ymax=1, color='silver', linewidth=0.8, linestyle='--')
axs.axhline(y=0, xmin=0, xmax=1, color='silver', linewidth=0.8, linestyle='--')
axs.axhline(y=1, xmin=0, xmax=1, color='silver', linewidth=0.8, linestyle='--')
# save
fig.savefig('./../fig/initial_state/MP_setting_NNPS_{}.png'.format(file_name), \
    format='png', dpi=300, transparent=False)
fig.savefig('./../fig/initial_state/MP_setting_NNPS_{}.pdf'.format(file_name), \
    format='pdf', transparent=True)
# close
plt.close()
print('[Message: 6/7] plot_MP_setting_NNPS has done.')

print('[Message: 7/7] plot_initial_state   has completed !')
end_time = time.perf_counter()
print('-----------------------------------------------------')
print('[Message] Total time : {:.2f} [s]'.format(end_time-start_time))
print('-----------------------------------------------------')
print('')


