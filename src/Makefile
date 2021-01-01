CC = gcc
CXX = g++
NVCC = nvcc

SRCS := $(wildcard */*.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRCS))

INCLUDES	= -I /usr/include -I ./src
LIBS		= -L/usr/lib/x86_64-linux-gnu -lboost_program_options -lpthread -lgeos
CPPFLAGS	= -g -std=c++14 -Wno-cpp

CUDALIBS	= -L/usr/local/cuda/lib64 -lcuda -lcudart

GEOMETRY_SRCS := $(wildcard geometry/*.cpp)
GEOMETRY_OBJS := $(patsubst %.cpp,%.o,$(GEOMETRY_SRCS))

GEOS_SRCS := $(wildcard geos/*.cpp)
GEOS_OBJS := $(patsubst %.cpp,%.o,$(GEOS_SRCS))

TRIANGULATE_SRCS := $(wildcard triangulate/*.cpp)
TRIANGULATE_OBJS := $(patsubst %.cpp,%.o,$(TRIANGULATE_SRCS))

#compile all the cu files
CU_SRCS := $(wildcard cuda/*.cu)
CU_OBJS := $(patsubst %.cu,%_cu.o,$(CU_SRCS))

#%_cu.o: %.cu
#	$(NVCC) -c $(CPPFLAGS) $(INCLUDES) -o $@ $<	
	
%.o:	%.cpp
	$(CXX) -c $(CPPFLAGS) $(INCLUDES) -o $@ $<	

all:	polygon_in_polygon point_polygon_distance point_in_polygon resque_polygon_in_polygon resque_point_polygon_distance resque_point_in_polygon \
		partition flatness pixeldensity distance lookup \
		load geninsert geninsert_point loadpoint


# for macro queries

polygon_in_polygon:	query/polygon_in_polygon.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
point_polygon_distance:	query/point_polygon_distance.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
point_in_polygon:	query/point_in_polygon.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	

geos_polygon_in_polygon:	query/geos_polygon_in_polygon.o $(GEOMETRY_OBJS) geos/GEOSTool.o
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
geos_point_polygon_distance:	query/geos_point_polygon_distance.o $(GEOMETRY_OBJS) geos/GEOSTool.o
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
geos_point_in_polygon:	query/geos_point_in_polygon.o $(GEOMETRY_OBJS) geos/GEOSTool.o
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
	
ms_point_in_polygon:	query/ms_point_in_polygon.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 

# utility

distance:	query/distance.o $(GEOMETRY_OBJS) 
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
lookup:	query/lookup.o $(GEOMETRY_OBJS) 
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
debug_distance:	query/debug_distance.o $(GEOMETRY_OBJS) 
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
debug_contain:	query/debug_contain.o $(GEOMETRY_OBJS) 
	$(CXX) -o ../build/$@ $^ $(LIBS) 


triangulate:	query/triangulate.o $(TRIANGULATE_OBJS) $(GEOMETRY_OBJS) 
	$(CXX) -o ../build/$@ $^ $(LIBS) 
# micro tests
partition:	stats/partition.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
flatness:	stats/flatness.o $(GEOMETRY_OBJS) 
	$(CXX) -o ../build/$@ $^ $(LIBS) 

pixeldensity:	stats/pixeldensity.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	


# for data extraction
load:	datagen/loadwkt.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 

loadpoint:	datagen/loadwkt_point.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
geninsert:	datagen/geninsert.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
geninsert_point:	datagen/geninsert_point.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 


	
clean:
	rm -fr ../build/* $(OBJS) $(CU_OBJS)
