#!gmake
##################################################
#
# Makefile for test for libmdlogp and libmdlogd
#
# $Id$
#
###################################################

include ../config.mk

DEVEL = ../../

AA_PATH = $(DEVEL)

INCLUDE = $(addprefix -I, $(AA_PATH)atomanalysis/include $(AA_PATH)atomanalysis/src ../inc $(DEVEL)libpka/inc )
INCLUDE += -I../inc -I../../grid_ng/inc -I../../grid_ng/libutil -I$(DEVEL)libbasstat/inc -I../../mini_mizer/lib_mmol/ \
	-I../../inchi/include/

CFLAGS += $(INCLUDE)
CXXFLAGS += $(INCLUDE)

ifeq ($(USEFLAP),yes)
  LIBRARY = -L../lib -llogdp -L ../../flap-db/lib/ -lflap -L$(AA_PATH)atomanalysis/lib -lAtomAnalysis \
	-L../../grid_ng/lib -lgrid -lgrin -lgrindutil -L../../libkekule/lib -lkekule
  LIBPKA = -L$(DEVEL)libpka/lib -lpka -L$(DEVEL)mini_mizer/lib -lmizer -lminimizer -lmmol -lmizer -lminimizer \
	 -L$(DEVEL)sqlite/.libs -lsqlite3 -L$(DEVEL)libbasstat -lBasStat -L$(DEVEL)mdtk/lib -lmdtk \
	 -L../../inchi/lib -linchi -lgfortran -lz -lstdc++ -lpthread
else
  LIBRARY = -L../lib -llogdp -L$(AA_PATH)atomanalysis/lib -lAtomAnalysis -lm \
	    -L../../grid_ng/lib -lgrid -lgrin -lgrindutil
  LIBPKA = -L$(DEVEL)libpka/lib -lpka -L$(DEVEL)mini_mizer/lib -lmizer -lminimizer -lmmol -lmizer -lminimizer \
           -L$(DEVEL)sqlite/.libs -lsqlite3 -L$(DEVEL)libbasstat -lBasStat -L$(DEVEL)mdtk/lib -lmdtk -lstdc++
endif


OBJ = plscoeffbld.o \
      common.o

all: plscoeffbld

plscoeffbld: $(OBJ)
	$(CC) -fopenmp -o $@ $(OBJ) $(LIBRARY) $(LIBPKA)
	cp -f $@ ../bin/

clean:
	-rm -f *.o plscoeffbld
