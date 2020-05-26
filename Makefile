SRCS := $(wildcard */*.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRCS))


INCFLAGS	= -I /usr/include -I ./src 
LIBS		= -L/usr/lib/x86_64-linux-gnu -lgeos -lspatialindex -lboost_program_options -lpthread
CPPFLAGS	= -g

all:	partitioner polygon

partitioner:	src/partitioner.o src/MyPolygon.o src/MyMultiPolygon.o src/Pixel.o
	$(CXX) -o build/$@ $^ $(LIBS) 

polygon:	src/Parser.o src/MyPolygon.o src/MyMultiPolygon.o src/Pixel.o
	$(CXX) -o build/$@ $^ $(LIBS) 
	
%.o:	%.cpp
	$(CXX) -c $(CFLAGS) $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

clean:
	rm -fr resque_2d $(OBJS)
