#!gmake
##################################################
#
# Makefile for test for libmdlogp and libmdlogd
#
# $Id$
#
###################################################

include ./config.mk

DEVEL = /home/redo/Sources/reposMD/

AA_PATH = $(DEVEL)

INCLUDE = $(addprefix -I, $(DEVEL)atomanalysis/include $(DEVEL)atomanalysis/src $(DEVEL)liblogd/inc $(DEVEL)libpka/inc )
INCLUDE += -I$(DEVEL)liblogd/inc -I$(DEVEL)grid_ng/inc -I$(DEVEL)grid_ng/libutil -I$(DEVEL)libbasstat/inc -I$(DEVEL)mini_mizer/lib_mmol/ \
	-I$(DEVEL)inchi/include/

CFLAGS += $(INCLUDE)
CXXFLAGS += $(INCLUDE)

ifeq ($(USEFLAP),yes)
  LIBRARY = -L../lib -llogdp -L $(DEVEL)flap-db/lib/ -lflap -L$(DEVEL)atomanalysis/lib -lAtomAnalysis \
	-L$(DEVEL)grid_ng/lib -lgrid -lgrin -lgrindutil -L$(DEVEL)libkekule/lib -lkekule
  LIBPKA = -L$(DEVEL)libpka/lib -lpka -L$(DEVEL)mini_mizer/lib -lmizer -lminimizer -lmmol -lmizer -lminimizer \
	 -L$(DEVEL)sqlite/ -lsqlite -L$(DEVEL)libbasstat -lBasStat -L$(DEVEL)mdtk/lib -lmdtk \
	 -L../../inchi/lib -linchi -lgfortran -lz -lstdc++ -lpthread
else
  LIBRARY = -L$(DEVEL)/liblogd/lib -llogdp -L$(DEVEL)atomanalysis/lib -lAtomAnalysis -lm \
	    -L$(DEVEL)grid_ng/lib -lgrid -lgrin -lgrindutil
  LIBPKA = -L$(DEVEL)libpka/lib -lpka -L$(DEVEL)mini_mizer/lib -lmizer -lminimizer -lmmol -lmizer -lminimizer \
           -L$(DEVEL)sqlite/ -lsqlite -L$(DEVEL)libbasstat -lBasStat -L$(DEVEL)mdtk/lib -lmdtk -lstdc++
endif


OBJ = plscoeffbld.o \
      common.o

OBJF = featuresbuild.o \
      common.o

all: plscoeffbld featuresbuild

plscoeffbld: $(OBJ)
	$(CC) -fopenmp -o $@ $(OBJ) $(LIBRARY) $(LIBPKA)

featuresbuild: $(OBJF)
	$(CC) -o $@ $(OBJF) $(LIBRARY) $(LIBPKA)

clean:
	-rm -f *.o plscoeffbld
