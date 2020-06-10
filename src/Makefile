CC = gcc
CXX = g++
NVCC = nvcc

SRCS := $(wildcard */*.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRCS))

INCLUDES	= -I /usr/include -I ./src
LIBS		= -L/usr/lib/x86_64-linux-gnu -lboost_program_options -lpthread
CPPFLAGS	= -g -std=c++14

CUDALIBS	= -L/usr/local/cuda/lib64 -lcuda -lcudart

GEOMETRY_SRCS := $(wildcard geometry/*.cpp)
GEOMETRY_OBJS := $(patsubst %.cpp,%.o,$(GEOMETRY_SRCS))

#compile all the cu files
CU_SRCS := $(wildcard cuda/*.cu)
CU_OBJS := $(patsubst %.cu,%_cu.o,$(CU_SRCS))

%_cu.o: %.cu
	$(NVCC) -c $(CPPFLAGS) $(INCLUDES) -o $@ $<	
%.o:	%.cpp
	$(CXX) -c $(CPPFLAGS) $(INCLUDES) -o $@ $<	

all:	partition query load geninsert

# for queries
partition:	query/partition.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 

polygon_in_polygon:	query/polygon_in_polygon.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
point_polygon_distance:	query/point_polygon_distance.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
point_in_polygon:	query/point_in_polygon.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
partstat:	query/partstat.o $(GEOMETRY_OBJS) $(CU_OBJS)
	$(CXX) -o ../build/$@ $^ $(CUDALIBS) $(LIBS) 

pixeldensity:	query/pixeldensity.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
distance:	query/distance.o $(GEOMETRY_OBJS) $(CU_OBJS)
	$(CXX) -o ../build/$@ $^ $(CUDALIBS) $(LIBS) 
	
lookup:	query/lookup.o $(GEOMETRY_OBJS) 
	$(CXX) -o ../build/$@ $^ $(LIBS) 

# for data extraction
load:	datagen/loadwkt.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
geninsert:	datagen/geninsert.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
geninsert_point:	datagen/geninsert_point.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 

loadpoint:	datagen/loadwkt_point.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
clean:
	rm -fr ../build/* $(OBJS) $(CU_OBJS)
