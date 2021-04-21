# Variables for directory paths
BDIR = ./build/
SDIR = ./src/streamLines/
LDIR = ./lib/
IDIR = ./lib/inc/

# Main target is defined in default
default: libVisual3.a

# Create libVisual3.a static library
libVisual3.a: $(LDIR)libVisual3.a StreamLines.o
	ar -r $(LDIR)libVisual3.a $(BDIR)StreamLines.o
	@echo "Library created successfully!!"

# Create modified StreamLines.o object file
StreamLines.o: $(SDIR)StreamLines.f
	gfortran -w -c -I $(IDIR) $(SDIR)StreamLines.f \
		 -o $(BDIR)StreamLines.o

# Clean target. Phony helps not to run it
.PHONY: clean

clean:
	rm -f $(BDIR)StreamLines.o
	@cp $(LDIR)libVisual3OLD.a $(LDIR)libVisual3.a
	@sudo sh -c "echo 1 > /proc/sys/vm/drop_caches"
	@echo "Library file is set to Visual3 Default"
	@echo "Cleared memory!!"
