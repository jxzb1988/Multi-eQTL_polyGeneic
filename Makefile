CC=g++
DIC=$(PWD)
CFLAGS=-c -Wall -g -O3  -I $(DIC)  
LDFLAGS= -DMKL_ILP64 -m64  -fopenmp -I $(DIC)/armadillo/include/ -DARMA_DONT_USE_WRAPPER -llapack -lblas -lgslcblas  -lgsl 
SOURCES1=MeQTLPolyG.cpp PostCal.cpp Util.cpp TopKSNP.cpp 
EXECUTABLE1=MeQTLPolyG
	
$(EXECUTABLE1): $(SOURCES1) 
	$(CC) $(SOURCES1)   $(LDFLAGS) -o $@
clean:
	rm MeQTLPolyG
