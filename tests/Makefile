SHELL = /bin/bash
# CXX = g++ 
CFLAGS = -Wextra -pedantic -Ofast -std=c++17 -Wall
INCLUDE = ../include


all: bbwt-test bbwt-console-test


bbwt-test: bbwt-test.cpp ${INCLUDE}/bbwt.hpp ${INCLUDE}/bbwt_internal.hpp ${INCLUDE}/lyndon.hpp ${INCLUDE}/BitVector.hpp
	${CXX} ${CFLAGS} -o bbwt-test bbwt-test.cpp -I${INCLUDE}
	
bbwt-console-test: bbwt-console-test.cpp ${INCLUDE}/bbwt.hpp ${INCLUDE}/bbwt_internal.hpp ${INCLUDE}/lyndon.hpp ${INCLUDE}/BitVector.hpp
	${CXX} ${CFLAGS} -o bbwt-console-test bbwt-console-test.cpp -I${INCLUDE}


clean:
	-${RM} bbwt-test bbwt-console-test
distclean: clean
	-${RM} bbwt-test bbwt-console-test

