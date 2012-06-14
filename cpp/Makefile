all: genefinder splitter

CXXFLAGS=-O2

genefinder: genefinder.o bindir
	c++ -o ../bin-cpp/genefinder genefinder.o
    
genefinder.o: genefinder.cpp
	c++ -c ${CXXFLAGS} genefinder.cpp
        
splitter: splitter.o bindir
	c++ -o ../bin-cpp/splitter splitter.o

splitter.o: splitter.cpp
	c++ -c ${CXXFLAGS} splitter.cpp

bindir:
	mkdir -p ../bin-cpp

clean:
	rm -rf *.o
    
distclean:
	rm -rf *.o
	rm -rf ../bin-cpp/genefinder
	rm -rf ../bin-cpp/splitter

    
