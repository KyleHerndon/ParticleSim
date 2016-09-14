#include <vector>
#include <cmath>
#include <stdlib.h>

#include "simulation.h"



int main() {
	int n = 2; // num_atoms = 8
	int dimensions = 3;
	double mass = 48.0;
	double temperature = 0.728;
	double boxLength = 1.0;
	Simulation sim = Simulation(dimensions, boxLength);
	sim.initCubic(n, temperature, mass);
}
