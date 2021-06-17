# Makefile 
#
.SUFFIXES:
.SUFFIXES: .o .f90

include ./user_build_config

OBJS =	ufsLandNoahMPReplaceRestart.o
	
all:	ufsLandNoahMPReplaceRestart.exe

.f90.o:
	$(COMPILERF90) -c $(F90FLAGS) $(FREESOURCE) $(NETCDFMOD) $(*).f90

ufsLandNoahMPReplaceRestart.exe: $(OBJS)
	$(COMPILERF90) -o $(@) $(F90FLAGS) $(FREESOURCE) $(NETCDFMOD) $(OBJS) $(NETCDFLIB)

clean:
	rm -f *.o *.mod *.exe


#
# Dependencies:
#
