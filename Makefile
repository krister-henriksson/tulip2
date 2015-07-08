

# #######################################################
# Defaults:
# #######################################################
prefix = $(HOME)
inc_libutils = $(prefix)/include/libutils
inc_spglib   = $(prefix)/include/spglib
lib_libutils = $(prefix)/lib
lib_spglib   = $(prefix)/lib
BINDIR       = $(prefix)/bin
# #######################################################


.PHONY: default dynamic dynamic2 static static2 dirs clean

INCDIR = $(inc_libutils) $(inc_spglib)
LIBDIR = $(lib_libutils) $(lib_spglib)
LIBDIRPATH = $(lib_libutils)":"$(lib_spglib)


# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# Settings:

CC      = g++
WARN    = -Wstrict-aliasing #-Wall -Wextra 
STD     = -ansi -pedantic -std=c++98 
DEBUG   = -g 
PROF    = -pg 
OPT     = -O2 -fstrict-aliasing #-funsafe-loop-optimizations -funroll-loops
OPENMP  = -fopenmp 

CFLAGS  = -c $(WARN) $(STD) $(DEBUG) $(OPT) $(PROF) $(OPENMP)

INC     = $(addprefix -I, $(INCDIR))
LIB     = $(addprefix -L, $(LIBDIR))

LDFLAGS        = $(LIB) -lutils  -lsymspg  -lm -lrt
LDFLAGS_STATIC = $(LDFLAGS) -static 

EXE_DYNAMIC := bin/tulip
EXE_STATIC  := bin/tulip_static
EXECUTABLES = $(EXE_DYNAMIC) $(EXE_STATIC)

SRC    := $(wildcard src/*.cpp)

SOURCES = $(SRC)
OBJECTS = $(SOURCES:src/%.cpp=obj/%.o)
DEPS    = obj/make.dep

REBUILDABLES = $(OBJECTS) $(EXECUTABLES)



default: dirs note $(OBJECTS)


install: default
	@echo ""
	@echo " *** Installing under " $(prefix)
	@echo " *** Installing executables in " $(BINDIR)
	@echo ""
	@echo ' *** If running bash, put these lines into ~/.bashrc :'
	@echo 'export LD_LIBRARY_PATH='$(LIBDIRPATH)':$${LD_LIBRARY_PATH}'
	@echo 'export LD_RUN_PATH='$(LIBDIRPATH)':$${LD_RUN_PATH}'
	@echo ""
	$(CC) $(OBJECTS) $(OPENMP) -o $(EXE_DYNAMIC)  $(LDFLAGS) $(PROF)
	# strip -s $(EXE_DYNAMIC)
	cp $(EXE_DYNAMIC) $(BINDIR)/
	$(CC) $(OBJECTS) $(OPENMP) -o $(EXE_STATIC)  $(LDFLAGS_STATIC) $(PROF)
	#strip -s $(EXE_STATIC)
	cp $(EXE_STATIC) $(BINDIR)/



# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# Pattern rules:

obj/%.o: src/%.cpp
	$(CC) $(INC) $(STD) $(WARN) $(DEBUG) $(OPT) $(PROF) $(OPENMP) -c  $< -o $@

$(DEPS): $(SOURCES)
	$(CC) -MM $(INC) $(SOURCES) | sed 's/\(.*\.o\)/obj\/\1/g' > $(DEPS)

include $(DEPS)

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# Phony rules:

dirs:
	- mkdir -p obj
	- mkdir -p bin
	- mkdir -p $(INCDIR)
	- mkdir -p $(LIBDIR)
	- mkdir -p $(BINDIR)

note:
	@echo ""
	@echo " *** Using these settings:"
	@echo " *** Include files for 'libutils'    : "$(inc_libutils)
	@echo " *** Include files for 'spglib'      : "$(inc_spglib)
	@echo " *** Library directory for 'libutils': "$(inc_libutils)
	@echo " *** Library directory for 'spglib'  : "$(inc_spglib)
	@echo ""


clean:
	-mkdir obj
	-rm -f obj/*


# $@: the target filename.
# $*: the target filename without the file extension.
# $<: the first prerequisite filename.
# $^: the filenames of all the prerequisites, separated by spaces, discard duplicates.
# $+: similar to $^, but includes duplicates.
# $?: the names of all prerequisites that are newer than the target, separated by spaces.


