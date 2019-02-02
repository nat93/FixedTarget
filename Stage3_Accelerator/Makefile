ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
EXTRALIBS += $(shell $(ROOTSYS)/bin/root-config --libs)
LDLIBS += $(shell $(ROOTSYS)/bin/root-config --libs)

CXX           = g++
CXXFLAGS      = -g -Wall -fPIC -Wno-deprecated

NGLIBS         = -Wl,--no-as-needed $(ROOTGLIBS)
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

all: accelerator

accelerator: accelerator.cc obj/mag.o obj/mag.o
	$(CXX) $(CXXFLAGS) -o accelerator $(OUTLIB)/*.o $(GLIBS) $<

obj/mag.o: src/MagClass.C src/MagClass.h 
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)mag.o $<

clean:
	rm -f accelerator
	rm -f obj/*.o
	rm -f *~ src/*~
