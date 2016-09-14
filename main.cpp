#include <vector>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
 

#include "simulation.h"

static void print_particle(std::ofstream &output, Particle a) {
	for (unsigned i = 0; i < a.getPosition().size(); i++) {
		output << " " << a.getPosition()[i];
	}
	output << "\n";
}

static void print_particles(std::ofstream &output, std::vector<Particle> particles) {
	for (unsigned i = 0; i < particles.size(); i++) {
		output << i;
		print_particle(output, particles[i]);
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

	std::ofstream output;
	output.open("trajectory.xyz");

	output << n*n*n << "\n"; //header

	print_particles(output, sim.getParticles());
	for (int i = 0; i < 1000; i++) {
		sim.updateParticles(.01);
		print_particles(output, sim.getParticles());
	}
	output.close();
	
	std::cout << sim.totalEnergy() << "\n";
}