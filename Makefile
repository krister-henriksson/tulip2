
.PHONY: dynamic dynamic2 static static2 dirs clean

include vars.mk


# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# Settings:

CC      = g++ 
WARN    = #-Wall -Wextra #-Wstrict-aliasing 
STD     = -ansi -pedantic -std=c++98 
DEBUG   = -g 
PROF    = -pg 
OPT     = -O3 #-fstrict-aliasing 
OPENMP  = -fopenmp 
CFLAGS  = -c $(WARN) $(STD) $(DEBUG) $(OPT) $(PROF) $(OPENMP)

INC     = -I$(INCDIR)

LIB            = -L$(LIBDIR)
LDFLAGS        = $(LIB) -lrt -lm -lutils
LDFLAGS_STATIC = $(LDFLAGS) -static 


EXECUTABLE = tulip



SRC    = compound.cpp \
	compound-getprop.cpp \
	elem-iacs.cpp \
	helpfuns.cpp \
	latcalc.cpp \
	mdsettings.cpp \
	mdsystem.cpp \
	mdsystem-pot-abop.cpp \
	mdsystem-pot-eam.cpp \
	param-pot.cpp \
	potclasses-pot-abop.cpp \
	potclasses-pot-eam.cpp \
	potclasses-reppot.cpp \
	potinfo.cpp \
	potinfo-readeam.cpp \
	potinfo-readinfo.cpp \
	potinfo-readreppot.cpp \
	potinfo-readspecs.cpp \
	propfun.cpp \
	report.cpp \
	specs-fit-prop-pot.cpp \
	$(EXECUTABLE).cpp

SOURCES = $(addprefix src/,$(SRC))
OBJECTS = $(SOURCES:src/%.cpp=obj/%.o)

#OBJECTS    = $(SOURCES:.cpp=.o)

REBUILDABLES = $(OBJECTS) $(EXECUTABLE)










# --------------------------------------------------------------------------
# --------------------------------------------------------------------------

dynamic: dirs dynamic2 $(EXECUTABLE)

dynamic2: $(OBJECTS)
	$(CC) $(OBJECTS) $(OPENMP) -o $(EXECUTABLE)  $(LDFLAGS) $(PROF)
	#strip -s $(EXECUTABLE)
	cp $(EXECUTABLE) $(BINDIR)/

# --------------------------------------------------------------------------

static: dirs static2 $(EXECUTABLE)

static2: $(OBJECTS)
	$(CC) $(OBJECTS) $(OPENMP) -o $(EXECUTABLE)  $(LDFLAGS_STATIC) $(PROF)
	#strip -s $(EXECUTABLE)
	cp $(EXECUTABLE) $(BINDIR)/


# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# Pattern rules:

obj/%.o: src/%.cpp
	$(CC) $(INC) $(STD) $(WARN) $(DEBUG) $(OPT) $(PROF) $(OPENMP) -c  $< -o $@



# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# Phony rules:

dirs:
	- mkdir obj
	- mkdir -p $(INCDIR)
	- mkdir -p $(LIBDIR)
	- mkdir -p $(BINDIR)

clean:
	-rm $(REBUILDABLES)


# $@: the target filename.
# $*: the target filename without the file extension.
# $<: the first prerequisite filename.
# $^: the filenames of all the prerequisites, separated by spaces, discard duplicates.
# $+: similar to $^, but includes duplicates.
# $?: the names of all prerequisites that are newer than the target, separated by spaces.


