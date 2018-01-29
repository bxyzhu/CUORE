CC = gcc
CXX = g++

ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

CXXFLAGS += -fPIC
CXXFLAGS += $(ROOTCFLAGS)
LDFLAGS  += $(ROOTLIBS)

all: run-wenqin MegaGPXEvolution libWenqinFitter.so libGPXFitter.so

MegaGPXEvolution: MegaGPXEvolution.cc libGPXFitter.so
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o MegaGPXEvolution MegaGPXEvolution.cc libGPXFitter.so

run-wenqin: run-wenqin.cc libWenqinFitter.so
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o run-wenqin run-wenqin.cc libWenqinFitter.so

libWenqinFitter.so: WenqinFitter.cc
	$(CXX) -shared $(CXXFLAGS) $(LDFLAGS) -lRooFitCore -lRooFit -lRooStats -o libWenqinFitter.so WenqinFitter.cc

libGPXFitter.so: GPXFitter.cc
	$(CXX) -shared $(CXXFLAGS) $(LDFLAGS) -lRooFitCore -lRooFit -lRooStats -o libGPXFitter.so GPXFitter.cc

clean:
	@echo "Removing all Wenqins...";
	rm -f run-wenqin MegaGPXEvolution libWenqinFitter.so libGPXFitter.so
