CC = gcc
CXX = g++

ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

CXXFLAGS += -fPIC
CXXFLAGS += $(ROOTCFLAGS)
LDFLAGS  += $(ROOTLIBS)

all: run-wenqin libWenqinFitter.so libWenqinSimulFitter.so

run-wenqin: run-wenqin.cc libWenqinFitter.so 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o run-wenqin run-wenqin.cc libWenqinFitter.so

libWenqinFitter.so: WenqinFitter.cc
	$(CXX) -shared $(CXXFLAGS) $(LDFLAGS) -lRooFitCore -lRooFit -lRooStats -o libWenqinFitter.so WenqinFitter.cc

libWenqinSimulFitter.so: WenqinSimulFitter.cc
	$(CXX) -shared $(CXXFLAGS) $(LDFLAGS) -lRooFitCore -lRooFit -lRooStats -o libWenqinSimulFitter.so WenqinSimulFitter.cc

clean:
	@echo "Cleaning ...";
	rm -f run-wenqin libWenqinFitter.so