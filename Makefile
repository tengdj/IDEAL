SRCS := $(wildcard */*.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRCS))


INCFLAGS	= -I /usr/include -I ./src 
LIBS		= -L/usr/lib/x86_64-linux-gnu -lgeos -lspatialindex -lboost_program_options -lpthread
CPPFLAGS	= -g

all:	resque_2d polygon

resque_2d:	src/resque_2d.o src/MyPolygon.o
	$(CXX) -o $@ $^ $(LIBS) 

polygon:	src/MyPolygon.o src/Parser.o
	$(CXX) -o $@ $^ $(LIBS) 
	
%.o:	%.cpp
	$(CXX) -c $(CFLAGS) $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

clean:
	rm -fr resque_2d $(OBJS)
