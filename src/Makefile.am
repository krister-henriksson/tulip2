
# http://www.jezuk.co.uk/cgi-bin/view/articles/autoconfiscation-part-two
# https://www.lrde.epita.fr/~adl/dl/autotools-handout-4.pdf


AUTOMAKE_OPTIONS = foreign


bin_PROGRAMS = tulip tulip_static

tulip_SOURCES = compound.cpp \
compoundfit.cpp \
compoundfit-getprop.cpp \
compoundfit.hpp \
compoundfit-list.cpp \
compoundfit-list.hpp \
compound.hpp \
elem-iacs.cpp \
elem-iacs.hpp \
errors.hpp \
get-comp-prop.cpp \
get-comp-prop.hpp \
helpfuns.cpp \
helpfuns.hpp \
lattice-simple.cpp \
lattice-simple.hpp \
mdsettings.cpp \
mdsettings.hpp \
mdsystem.cpp \
mdsystem-force.cpp \
mdsystem.hpp \
mdsystem-pot-abop.cpp \
mdsystem-pot-eam.cpp \
mdsystem-relax.cpp \
param-pot.cpp \
param-pot.hpp \
physconst.hpp \
potclasses.hpp \
potclasses-pot-abop.cpp \
potclasses-pot-eam.cpp \
potclasses-reppot.cpp \
potinfo.cpp \
potinfo.hpp \
potinfo-readeam.cpp \
potinfo-readinfo.cpp \
potinfo-readreppot.cpp \
potinfo-readspecs.cpp \
propfun.cpp \
propfun.hpp \
report.cpp \
report.hpp \
specs-fit-prop-pot.cpp \
specs-fit-prop-pot.hpp \
tulip.cpp

tulip_static_SOURCES = $(tulip_SOURCES)


ma_inc  = -I$(srcdir)
ma_warn =  -Wstrict-aliasing # -Wall -Wextra
ma_std  = -ansi -pedantic -std=c++98
ma_opt  = -O2
ma_dbg  = -g
ma_omp  = -fopenmp


tulip_CPPFLAGS = $(ma_inc) $(ma_warn) $(ma_std) $(ma_opt) $(ma_dbg) $(ma_omp)
tulip_static_CPPFLAGS = $(ma_inc) $(ma_warn) $(ma_std) $(ma_opt) $(ma_dbg) $(ma_omp)

tulip_LDADD = -lutils -lsymspg -lm -lrt
tulip_LDFLAGS = $(ma_omp)

tulip_static_LDADD = -lutils -lsymspg -lm -lrt
tulip_static_LDFLAGS = $(ma_omp) -static

