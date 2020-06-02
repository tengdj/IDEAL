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

partition:	tool/partition.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 

query:	tool/query.o $(GEOMETRY_OBJS) $(CU_OBJS)
	$(CXX) -o ../build/$@ $^ $(CUDALIBS) $(LIBS) 
	
load:	tool/loadwkt.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
geninsert:	tool/geninsert.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 

partperf:	tool/partperf.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 

clean:
	rm -fr ../build/* $(OBJS) $(CU_OBJS)
