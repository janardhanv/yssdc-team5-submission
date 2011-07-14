all: rnachecker

CXXFLAGS=-O2

rnachecker: checking.o main.o gbk_reader.o answer_reader.o
	mkdir -p ../bin-cpp
	c++ -o ../bin-cpp/rnachecker checking.o main.o gbk_reader.o answer_reader.o
    
checking.o: checking.cpp checking.h
	c++ -c ${CXXFLAGS} checking.cpp
    
gbk_reader.o: gbk_reader.cpp gbk_reader.h
	c++ -c ${CXXFLAGS} gbk_reader.cpp

answer_reader.o: answer_reader.cpp answer_reader.h
	c++ -c ${CXXFLAGS} answer_reader.cpp

main.o: main.cpp answer_reader.h gbk_reader.h
	c++ -c ${CXXFLAGS} main.cpp
    
clean:
	rm -rf *.o
    
distclean:
	rm -rf *.o
	rm -rf rnachecker

    
