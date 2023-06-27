# name of TARGET file
TARGET = start_calculation
# compiler
F90 = gfortran
# compile options
FFLAGS += -O -fopenmp -fconvert=big-endian -fimplicit-none -fbounds-check
#FFLAGS += -O -fconvert=big-endian -fimplicit-none -fbounds-check

# OBJESTS files
OBJECTS = \
	input.o \
	glb_setting.o \
	glb_arg.o \
	cal_EOS.o \
	initial_setting.o \
	NNPS.o \
	cal_kernel.o \
	set_time_step.o \
	cal_wall.o \
	cal_main.o \
	PS_scheme.o \
	cal_mapping.o \
	analysis.o \
	write_out.o \
	main.o

# MODULE files
COMMAN_MOD = \
	input.f90 \
	glb_setting.f90 \
	glb_arg.f90 \
	cal_EOS.f90 \
	initial_setting.f90 \
	NNPS.f90 \
	cal_kernel.f90 \
	set_time_step.f90 \
	cal_wall.f90 \
	cal_main.f90 \
	PS_scheme.f90 \
	cal_mapping.f90 \
	analysis.f90 \
	write_out.f90 \
	main.f90

#============= DO NOT CHANGE ==============#
# implicit none
.SUFFIXES :
# f90 -> o
.SUFFIXES : .o .f90
# compile each files
.f90.o:
	$(F90) -c $(FFLAGS) $<
# make binary file
$(TARGET) : $(OBJECTS)
	$(F90) $(FFLAGS) -o $@ $(OBJECTS)
	@echo ------
	@echo [Message] Makefile has been made ! 
	@echo [Message] put 'start_calculation'
# if MODULE files changed
$(OBJECTS) : $(COMMAN_MOD)

# make clean
.PHONY : clean
clean :
#	@rm -rf $(TARGET).exe $(OBJECTS) *.mod
	@del $(TARGET).exe $(OBJECTS) *.mod
	@echo [Message] $(TARGET).exe, *.o, and *.mod have been removed.