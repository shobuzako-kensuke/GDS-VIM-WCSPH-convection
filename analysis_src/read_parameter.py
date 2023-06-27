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
# read parameter
#=====================
# open and read parameter set
path_para = glob.glob('./../output_setting/for_python_*.dat')
data_para = open(path_para[0], 'r', encoding='UTF-8')
read_para = data_para.readlines()
data_para.close()
# remove '/n'
for i in range(len(read_para)):
    read_para[i] = read_para[i].strip()
# read each parameter
total_step       = int(float(read_para[0]))
num              = int(float(read_para[1]))
num_system       = int(float(read_para[2]))
total_num        = int(float(read_para[3]))
total_num_system = int(float(read_para[4]))
num_wall         = int(float(read_para[5]))
num_MP_side      = int(float(read_para[6]))
num_MP           = int(float(read_para[7]))
wall_thickness   = int(float(read_para[8]))
write_step       = int(float(read_para[9]))
RSST_VIM         = int(float(read_para[10]))
kernel_switch    = int(float(read_para[11]))
 
length         = float(read_para[12])
delta_tem      = float(read_para[13])
coe_smooth     = float(read_para[14])
delta_x        = float(read_para[15])
smooth_len     = float(read_para[16])
support_domain = float(read_para[17])
time_step      = float(read_para[18])
 
rho_0   = float(read_para[19])
alpha   = float(read_para[20])
eta_0   = float(read_para[21])
T_0     = float(read_para[22])
T_top   = float(read_para[23])
kappa   = float(read_para[24])
 
V_fin    = float(read_para[25])
V_fin_tb = float(read_para[26])
 
Ra   = float(read_para[27])
Pr   = float(read_para[28])
M_th = float(read_para[29])
SS   = float(read_para[30])
RSS  = float(read_para[31])