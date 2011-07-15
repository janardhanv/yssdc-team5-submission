all: genefinder

CXXFLAGS=-O2

genefinder: genefinder.o
	mkdir -p ../bin-cpp
	c++ -o ../bin-cpp/genefinder genefinder.o
    
genefinder.o: genefinder.cpp
	c++ -c ${CXXFLAGS} genefinder.cpp
        
clean:
	rm -rf *.o
    
distclean:
	rm -rf *.o
	rm -rf genefinder

    
