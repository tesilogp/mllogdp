CC = gcc
CXX = g++
CFLAGS = -Wall -W -O0 -g -I../inc -I../../atomanalysis/include \
	 -I../../atomanalysis/src -D_AA_TOOLS -I../../libpka/inc 
CXXFLAGS = -Wall -W -O0 -g -I../inc -I../../atomanalysis/include \
	 -I../../atomanalysis/src -D_AA_TOOLS -I../../libpka/inc \
	 -I../../mdtk

# force to use the libpka new api even if LIBPKA_VERSION is not defined
FORCENEWCLPKAAPI = no

# use flap extenision 
USEFLAP = no

ifeq ($(FORCENEWCLPKAAPI),yes)
CFLAGS += -DFORCE_NEW_C_LIPBPKA_API  
endif

MAJOR = 0
MINOR = 1
PATCH = 0

ifeq ($(USEFLAP),yes)
  CXXFLAGS += -DMAJOR=$(MAJOR) -DMINOR=$(MINOR) -DPATCH=$(PATCH) -DUSE_ALSO_FLAP
  CFLAGS += -DMAJOR=$(MAJOR) -DMINOR=$(MINOR) -DPATCH=$(PATCH) -DUSE_ALSO_FLAP
else
  CXXFLAGS += -DMAJOR=$(MAJOR) -DMINOR=$(MINOR) -DPATCH=$(PATCH) 
  CFLAGS += -DMAJOR=$(MAJOR) -DMINOR=$(MINOR) -DPATCH=$(PATCH) 
endif

