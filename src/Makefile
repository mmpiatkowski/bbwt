SHELL = /bin/bash
# CXX = g++
CFLAGS = -Wextra -pedantic -Ofast -std=c++17 -Wall
INCLUDE = ../include


all: bbwt bbwt-console csa-console


bbwt: bbwt-main.cpp ${INCLUDE}/bbwt.hpp ${INCLUDE}/bbwt_internal.hpp ${INCLUDE}/lyndon.hpp ${INCLUDE}/BitVector.hpp
	${CXX} ${CFLAGS} -o bbwt bbwt-main.cpp -I${INCLUDE}
	
bbwt-console: bbwt-console.cpp ${INCLUDE}/bbwt.hpp ${INCLUDE}/bbwt_internal.hpp ${INCLUDE}/lyndon.hpp ${INCLUDE}/BitVector.hpp
	${CXX} ${CFLAGS} -o bbwt-console bbwt-console.cpp -I${INCLUDE}

csa-console: csa-console.cpp ${INCLUDE}/bbwt.hpp ${INCLUDE}/bbwt_internal.hpp ${INCLUDE}/lyndon.hpp ${INCLUDE}/BitVector.hpp
	${CXX} ${CFLAGS} -o csa-console csa-console.cpp -I${INCLUDE}

clean:
	-${RM} bbwt bbwt-console csa-console
distclean: clean
	-${RM} bbwt bbwt-console csa-console

