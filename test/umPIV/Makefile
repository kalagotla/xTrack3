#
LIB = ../../lib

all: init gridData.o flowData.o

init:	init.o FORCE
	gfortran -o init init.o\
                 -L$(LIB) -lVisual3 -lGLU -lGL \
		 -L/usr/X11R6/lib64 -lX11 -lXext

gridData.o: gridData.f90
	@echo "Creating Grid File"
	@gfortran gridData.f90
	@./a.out
	@echo "Grid File... Done!"

flowData.o: flowData.f90
	@echo "Creating Flow File"
	@gfortran flowData.f90
	@./a.out
	@echo "Flow File... Done!"

init.o:	init.f
	gfortran -g -c -O -fno-automatic init.f

# Force init to link to lib whenever make is run
FORCE:

.PHONY: clean

clean:
	@rm -f init init.o gridData.o flowData.o a.out
	@echo "Removed object files"
