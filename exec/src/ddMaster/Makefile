INCFLAGS = -I/usr/local/include/ -I../

CPP = g++
CPPFLAGS = -g -O3 $(INCFLAGS)  -fopenmp -Wno-strict-aliasing
LINKERFLAGS = -lz
DEBUGFLAGS = -g -ggdb $(INCFLAGS)
HEADERS=$(wildcard *.h**)

ddMaster: ddMaster.cpp $(HEADERS)
	$(CPP) $(CPPFLAGS) ddMaster.cpp -o ddMaster $(LINKERFLAGS)

clean:
	@rm -rf ddMaster
