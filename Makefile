CC = g++
CFLAGS = -Wall -std=c++11

all: simulation

simulation: particle.o simulation.o main.o
	$(CC) $(CFLAGS) -o simulation particle.o simulation.o main.o

main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp

particle.o: particle.cpp particle.h
	$(CC) $(CFLAGS) -c particle.cpp

simulation.o: simulation.cpp simulation.h
<<<<<<< HEAD
	$(CC) $(CFLAGS) -c simulation.cpp
=======
	g++ -Wall -c simulation.cpp -std=c++11

clean:
	rm particle.o simulation.o
>>>>>>> 5ec2bc1432ea3a393875fdc119c74c41ffab69b1
