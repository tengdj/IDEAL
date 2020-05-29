SRCS := $(wildcard */*.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRCS))


LIBS		= -L/usr/lib/x86_64-linux-gnu -lboost_program_options -lpthread
CPPFLAGS	= -g -I /usr/include -I ./src 

GEOMETRY_SRCS := $(wildcard geometry/*.cpp)
GEOMETRY_OBJS := $(patsubst %.cpp,%.o,$(GEOMETRY_SRCS))

all:	partition query load geninsert

partition:	tool/partition.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 

query:	tool/query.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
load:	tool/loadwkt.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
geninsert:	tool/geninsert.o $(GEOMETRY_OBJS)
	$(CXX) -o ../build/$@ $^ $(LIBS) 
	
%.o:	%.cpp
	$(CXX) -c $(CFLAGS) $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

teng:
	echo $(GEOMETRY_SRCS)
clean:
	rm -fr ../build/* $(OBJS)
