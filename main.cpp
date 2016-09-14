#include <vector>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <math.h>
 

#include "simulation.h"

static void print_particle(Particle a) {
	std::cout << "(";
	for (double i : a.getPosition()) {
		std::cout << i << " ";
	}
	std::cout << ")\n";
}

static void print_particles(std::vector<Particle> particles) {
	for (Particle &particle : particles) {
		print_particle(particle);
	}
}

int main() {
	int n = 4; // num_atoms = 64
	int dimensions = 3;
	double mass = 48.0;
	double temperature = 0.728;
	double boxLength = 4.2323167;
	Simulation sim = Simulation(dimensions, boxLength);
	sim.initCubic(n, temperature, mass);
	std::cout << sim.totalEnergy() << "\n";
	print_particles(std::vector<Particle> ());
	//print_particles(sim.getParticles());
	for (int i = 0; i < 1000; i++) {
		sim.updateParticles(0.001);
	}
	std::cout << sim.totalEnergy() << "\n";
	//print_particles(sim.getParticles());
}