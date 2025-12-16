CXX = clang++
CXXFLAGS = -g -O3 -Wall -fopenmp=libomp
LDFLAGS = -lm -lhdf5 -lhdf5_cpp -lstdc++

all: beam.o mesh.o ray_trace.o main.o cbet.o
	$(CXX) $(CXXFLAGS) -o cbet beam.o mesh.o ray_trace.o main.o cbet.o $(LDFLAGS)

cbet.o: cbet.cpp beam.hpp mesh.hpp consts.hpp
	$(CXX) $(CXXFLAGS) -c cbet.cpp $(LDFLAGS)

main.o: main.cpp beam.hpp mesh.hpp cbet.hpp consts.hpp ray_trace.hpp
	$(CXX) $(CXXFLAGS) -c main.cpp $(LDFLAGS)

beam.o: beam.cpp beam.hpp consts.hpp
	$(CXX) $(CXXFLAGS) -c beam.cpp $(LDFLAGS)

mesh.o: mesh.cpp mesh.hpp consts.hpp
	$(CXX) $(CXXFLAGS) -c mesh.cpp $(LDFLAGS)

ray_trace.o: ray_trace.cpp ray_trace.hpp consts.hpp beam.hpp mesh.hpp
	$(CXX) $(CXXFLAGS) -c ray_trace.cpp $(LDFLAGS)
