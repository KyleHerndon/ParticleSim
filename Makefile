all: particle.o simulation.o

particle.o: particle.cpp particle.h
	g++ -Wall -c particle.cpp -std=c++11

simulation.o: simulation.cpp simulation.h
	g++ -Wall -c simulation.cpp -std=c++11

clean:
	rm particle.o simulation.o
