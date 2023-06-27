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

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#                           INPUT below                               #
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# if you want to delete files, input 1
output_setting_del     = 1         # files in output_setting
save_file_del          = 1         # files in save_file
fig_initial_del        = 1         # figures in initial_setting
fig_analysis_del       = 1         # figures in analysis
fig_movie_del          = 1         # figures in movie_material
fig_movie_wall_del     = 1         # figures in movie_material_wall
fig_movie_all_del      = 1         # figures in movie_material_all
fig_movie_all_wall_del = 1         # figures in movie_material_all_wall
fig_movie_PS_shift     = 1         # figures in movie_material_PS_shift
makefile_del           = 1         # makefile
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#                           INPUT end                                 #
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#---------------------
# NOT CHANGE
#---------------------
def remove_data(path_name):
    # get file name
    for i in path_name:
        # if file exists, remove_exe
        if os.path.isfile(i):
            os.remove(i)

#---------------------
# search files
#---------------------
if (output_setting_del == 1):
    path_name = glob.glob('./output_setting/*.dat')
    remove_data(path_name)
    print('[Message] Files   in output_setting          have been removed.')
if (save_file_del == 1):
    path_name_1 = glob.glob('./save_file/**/*.dat')
    path_name_2 = glob.glob('./save_file/*.dat')
    path_name = path_name_1 + path_name_2
    remove_data(path_name)
    print('[Message] Files   in save_file               have been removed.')
if (fig_initial_del == 1):
    path_name_1 = glob.glob('./fig/initial_state/*.pdf')
    path_name_2 = glob.glob('./fig/initial_state/*.png')
    path_name = path_name_1 + path_name_2
    remove_data(path_name)
    print('[Message] Figures in initial_setting         have been removed.')
if (fig_analysis_del == 1):
    path_name_1 = glob.glob('./fig/analysis/*.pdf')
    path_name_2 = glob.glob('./fig/analysis/*.png')
    path_name = path_name_1 + path_name_2
    remove_data(path_name)
    print('[Message] Figures in analysis                have been removed.')
if (fig_movie_del == 1):
    path_name = glob.glob('./fig/movie_material/*.png')
    remove_data(path_name)
    print('[Message] Figures in movie_material          have been removed.')
if (fig_movie_wall_del == 1):
    path_name = glob.glob('./fig/movie_material_wall/*.png')
    remove_data(path_name)
    print('[Message] Figures in movie_material_wall     have been removed.')
if (fig_movie_all_del == 1):
    path_name = glob.glob('./fig/movie_material_all/*.png')
    remove_data(path_name)
    print('[Message] Figures in movie_material_all      have been removed.')
if (fig_movie_all_wall_del == 1):
    path_name = glob.glob('./fig/movie_material_all_wall/*.png')
    remove_data(path_name)
    print('[Message] Figures in movie_material_all_wall have been removed.')
if (fig_movie_PS_shift == 1):
    path_name = glob.glob('./fig/PS_shift/*.png')
    remove_data(path_name)
    print('[Message] Figures in movie_material_PS_shift have been removed.')
if (makefile_del == 1):
    path_name_1 = glob.glob('./program_src/*.mod')
    path_name_2 = glob.glob('./program_src/*.o')
    path_name = path_name_1 + path_name_2
    remove_data(path_name)
    print('[Message] .mod .o in program_src             have been removed.')