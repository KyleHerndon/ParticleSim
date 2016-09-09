all: particle.o

particle.o: particle.cpp particle.h
	g++ -Wall -c particle.cpp
