fext:=f08
FC:=gfortran

FFLAGS:=
main_program:=main.$(fext)
sourcefiles:=$(wildcard *.$(fext))	#all files ending in specifiead fortran extension
objectfiles:=$(patsubst %.$(fext),%.o,$(sourcefiles))	#their corresponding module names
modules:=$(patsubst %.$(fext),%.mod,$(sourcefiles))	#their corresponding module names
modules:=$(filter-out $(main_program),$(modules))	#filter out name of main-file since it produces no .mod file
#
build:$(objectfiles)
	$(FC) -o main $(objectfiles)
	
%.o:%.$(fext)
	$(FC) $(FFLAGS) -c $<
