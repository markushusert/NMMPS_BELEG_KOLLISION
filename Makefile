fext:=f08
FC:=gfortran
make_depend_file:=make_depend
FFLAGS:=-g
main_program:=main.$(fext)
sourcefiles:=$(wildcard *.$(fext))	#all files ending in specifiead fortran extension
objectfiles:=$(patsubst %.$(fext),%.o,$(sourcefiles))	#their corresponding module names
modules:=$(patsubst %.$(fext),%.mod,$(sourcefiles))	#their corresponding module names
modules:=$(filter-out $(main_program),$(modules))	#filter out name of main-file since it produces no .mod file
#

build:$(objectfiles) $(modules)
	$(FC) -o main $(objectfiles)
	@echo -----------succesfully compiled program-----------
-include $(make_depend_file)
%.o %.mod:%.$(fext)
	$(FC) $(FFLAGS) -c $<
	#echo building target $@, prerequisites $^
	@touch $@

clear:
	rm -f $(objectfiles) $(modules) $(make_depend_file)
depend:$(sourcefiles)
	rm -f $(make_depend_file)
	makedepf90 $(sourcefiles) > $(make_depend_file)
	
