CC = gcc
CXX = g++
NVCC = nvcc

SRCS := $(wildcard */*.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRCS))

INCLUDES	= -I /usr/include -I include -I /usr/local/include
#-lgeos_c  
LIBS		= -L/usr/lib/x86_64-linux-gnu -L/usr/local/lib -lboost_program_options -lpthread -fopenmp -lgeos  
CPPFLAGS	= -g -std=c++14 -Wno-cpp -fopenmp
NVCCFLAGS	= -g -std=c++14 -DUSE_GPU

GEOMETRY_SRCS := $(wildcard geometry/*.cpp)
GEOMETRY_SRCS += $(wildcard geometry/triangulate/*.cpp)
GEOMETRY_SRCS += $(wildcard storage/*.cpp)
GEOMETRY_OBJS := $(patsubst %.cpp,%.o,$(GEOMETRY_SRCS))

PARTITION_SRCS := $(wildcard partition/*.cpp)
PARTITION_OBJS := $(patsubst %.cpp,%.o,$(PARTITION_SRCS))

#compile all the cu files
CU_SRCS := $(wildcard cuda/*.cu)
CU_OBJS := $(patsubst %.cu,%_cu.o,$(CU_SRCS))

ifdef USE_GPU
	CPPFLAGS += -DUSE_GPU
	LIBS  += -L/usr/local/cuda/lib64 -lcuda -lcudart
	GEOMETRY_OBJS += ${CU_OBJS}
endif


%_cu.o: %.cu
	$(NVCC) -c $(NVCCFLAGS) $(INCLUDES) -o $@ $<	
	
%.o:	%.cpp
	$(CXX) -c $(CPPFLAGS) $(INCLUDES) -o $@ $<	

all:	contain_polygon within_polygon


# for macro queries

contain_polygon:	query/contain_polygon.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
within_polygon:	query/within_polygon.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
within:	query/within.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
contain:	query/contain.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
print_raster:	test/print_raster.o $(GEOMETRY_OBJS)
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

test:	test/test.o $(GEOMETRY_OBJS)
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

loadpais:	datagen/loadwkt_pais.o $(GEOMETRY_OBJS) 
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

print:	test/print.o $(GEOMETRY_OBJS) $(PARTITION_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
printspace:	test/print_spaces.o $(GEOMETRY_OBJS) $(PARTITION_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
sample:	test/sample.o $(GEOMETRY_OBJS) $(PARTITION_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 

subset:	test/subset.o $(GEOMETRY_OBJS) $(PARTITION_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
partition_stats:	test/partition_stats.o $(GEOMETRY_OBJS) $(PARTITION_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 

clean:
	rm -fr ../build/* $(OBJS) $(CU_OBJS)
