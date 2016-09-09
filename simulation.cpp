#include "simulation.h"

#include <cmath>

Simulation::Simulation();

std::vector<double> getDisplacement(Particle a, Particle b) {
	std::vector<double> v = std::vector<double>(dim, 0);
	for (int i = 0; i < dim; i++) {
		double temp = a.getPosition()[i] - b.getPosition()[i];
		if (temp < (boxLength / 2) && temp > -(boxLength / 2)) {
			v.assign(i, temp)
		} else {
			v.assign(i, (temp > 0) ? boxLength-temp : -boxLength+temp);
		}
	}
}
double getDistance(Particle a, Particle b) {
	double ret = 0;
	for (double component : getDisplacement(a, b)) {
		ret += std::pow(component, 2);
	}
	return std::sqrt(ret);
}
std::vector<double> positionInBox(Particle a);