
.PHONY: dynamic dynamic2 static static2 dirs clean

include vars.mk


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

LIB            = $(addprefix -L, $(LIBDIR))
LDFLAGS        = $(LIB) -lutils  -lsymspg  -lm -lrt
#LDFLAGS        = $(LIB) -lutils -lm -lrt
LDFLAGS_STATIC = $(LDFLAGS) -static 


EXECUTABLE := bin/tulip


SRC    := $(wildcard src/*.cpp)

SOURCES = $(SRC)
OBJECTS = $(SOURCES:src/%.cpp=obj/%.o)
DEPS = obj/make.dep

REBUILDABLES = $(OBJECTS) $(EXECUTABLE)




# --------------------------------------------------------------------------
# --------------------------------------------------------------------------

dynamic: dirs dynamic2 $(EXECUTABLE)

dynamic2: $(OBJECTS)
	$(CC) $(OBJECTS) $(OPENMP) -o $(EXECUTABLE)  $(LDFLAGS) $(PROF)
	# strip -s $(EXECUTABLE)
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

clean:
	- rm -f $(REBUILDABLES) $(DEPS)


# $@: the target filename.
# $*: the target filename without the file extension.
# $<: the first prerequisite filename.
# $^: the filenames of all the prerequisites, separated by spaces, discard duplicates.
# $+: similar to $^, but includes duplicates.
# $?: the names of all prerequisites that are newer than the target, separated by spaces.


