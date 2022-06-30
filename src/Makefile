CC = gcc
CXX = g++
NVCC = nvcc

SRCS := $(wildcard */*.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRCS))

INCLUDES	= -I /usr/include -I include
LIBS		= -L/usr/lib/x86_64-linux-gnu -lboost_program_options -lpthread -lgeos -fopenmp
CPPFLAGS	= -g -std=c++14 -Wno-cpp 

CUDALIBS	= -L/usr/local/cuda/lib64 -lcuda -lcudart

GEOMETRY_SRCS := $(wildcard geometry/*.cpp)
GEOMETRY_SRCS += $(wildcard geometry/triangulate/*.cpp)
GEOMETRY_SRCS += $(wildcard storage/*.cpp)
GEOMETRY_OBJS := $(patsubst %.cpp,%.o,$(GEOMETRY_SRCS))

PARTITION_SRCS := $(wildcard partition/*.cpp)
PARTITION_OBJS := $(patsubst %.cpp,%.o,$(PARTITION_SRCS))

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

contain_polygon:	query/contain_polygon.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
within_polygon:	query/within_polygon.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
within:	query/within.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
contain:	query/contain.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	

geos_polygon_in_polygon:	query/geos_polygon_in_polygon.o $(GEOMETRY_OBJS) 
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
geos_point_polygon_distance:	query/geos_point_polygon_distance.o $(GEOMETRY_OBJS) 
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
geos_point_in_polygon:	query/geos_point_in_polygon.o $(GEOMETRY_OBJS) 
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
geos_polygon_polygon_within:	query/geos_polygon_polygon_within.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 

# utility

distance:	test/distance.o $(GEOMETRY_OBJS) 
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
lookup:	test/lookup.o $(GEOMETRY_OBJS) 
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
debug_distance:	test/debug_distance.o $(GEOMETRY_OBJS) 
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
debug_contain:	test/debug_contain.o $(GEOMETRY_OBJS) 
	$(CXX) -o ../build/$@ $^ $(LIBS) 

test:	test/within_polygon.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 	
	
# micro tests
#partition:	stats/partition.o $(GEOMETRY_OBJS) $(TRIANGULATE_OBJS) 
#	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
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

# for data partitioning

partition_analysis:	test/partition_analysis.o $(GEOMETRY_OBJS) $(PARTITION_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
genschema:	test/generate_partition.o $(GEOMETRY_OBJS) $(PARTITION_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 

printpoints:	test/print_points.o $(GEOMETRY_OBJS) $(PARTITION_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
sample:	test/sample.o $(GEOMETRY_OBJS) $(PARTITION_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 

subset:	test/subset.o $(GEOMETRY_OBJS) $(PARTITION_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
clean:
	rm -fr ../build/* $(OBJS) $(CU_OBJS)
