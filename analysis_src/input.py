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
import subprocess
rc('text', usetex=True)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
#                           INPUT below                               #
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#

# if exe, put 1
do_plot_movie_material = 1
do_plot_analysis       = 1
do_plot_initial_state  = 1

# input save file number you want to analyse
# see in ./../save_file/V_rms_save*.dat
# [note] put (the number - 1)
num_save_file = 50

#======================
# for movie materials
#======================
# movie switch (if exe, put 1)
plot_without_wall     = 1
plot_with_wall        = 0
plot_all_without_wall = 0
plot_all_with_wall    = 0
plot_PS_shift         = 0

# select physical quantity scattered in movie_materials
# (2,3,4,5,6,7,8,9) = (u,v,rho_S, pre, tem, eta, rho_sm, rho_G)
plot_phy = 6

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
#                           INPUT end                                 #
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#

if (do_plot_movie_material == 1):
    import plot_movie_material

if (do_plot_analysis == 1):  
    import plot_analysis

if (do_plot_initial_state == 1):
    import plot_initial_state    
