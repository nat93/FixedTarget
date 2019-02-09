ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
EXTRALIBS += $(shell $(ROOTSYS)/bin/root-config --libs)
LDLIBS += $(shell $(ROOTSYS)/bin/root-config --libs)

CXX           = g++
CXXFLAGS      = -g -Wall -fPIC -Wno-deprecated

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit

CXXFLAGS      += $(ROOTCFLAGS)
CXX           += -I./
LIBS           = $(ROOTLIBS) 

GLIBS          = $(filter-out -lNew, $(NGLIBS))

CXX	      += -I./obj/
OUTLIB	      = ./obj/
.SUFFIXES: .cc,.C
.PREFIXES: ./obj/

#----------------------------------------------------#

all: hardccbar1 hardccbar2

hardccbar1: hardccbar1.cc
	g++ $(CXXFLAGS) -o hardccbar1 hardccbar1.cc -I/home/andrii/pythia/pythia8185/include -L/home/andrii/pythia/pythia8185/lib/archive -lpythia8 -llhapdfdummy `root-config --cflags --glibs`

hardccbar2: hardccbar2.cc
	g++ $(CXXFLAGS) -o hardccbar2 hardccbar2.cc -I/home/andrii/pythia/pythia8185/include -L/home/andrii/pythia/pythia8185/lib/archive -lpythia8 -llhapdfdummy `root-config --cflags --glibs`


clean:
	rm -f *~
	rm -f src/*~
	rm -f $(OUTLIB)*.o
	rm -f hardccbar1
	rm -f hardccbar2
