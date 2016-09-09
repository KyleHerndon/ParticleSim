#include "simulation.h"

#include <cmath>

Simulation::Simulation(int d, double l) {
	dim = d;
	boxLength = l;
}

std::vector<double> Simulation::displacement(Particle a, Particle b) {
	std::vector<double> v = std::vector<double>(dim, 0);
	for (int i = 0; i < dim; i++) {
		double temp = a.getPosition()[i] - b.getPosition()[i];
		if (temp < (boxLength / 2) && temp > -(boxLength / 2)) {
			v[i] = temp;
		} else {
			v[i] = (temp > 0) ? boxLength-temp : -boxLength+temp;
		}
	}
	return v;
}

double Simulation::distance(Particle a, Particle b) {
	double ret = 0;
	for (double component : displacement(a, b)) {
		ret += std::pow(component, 2);
	}
	return std::sqrt(ret);
}

std::vector<double> Simulation::positionInBox(Particle a) {
	std::vector<double> newPos = std::vector<double>();
	for (double component : a.getPosition()) {
		newPos.push_back((component > (boxLength / 2) ? boxLength-component : component));
	}
	return newPos;
}

std::vector<double> Simulation::internalForce(Particle a) {
	std::vector<double> force = std::vector<double>(dim, 0);
	for (Particle particle : particles) {
		if (&particle != &a) {
			std::vector<double> temp = Simulation::force(a, particle);
			for (int i = 0; i < dim; i++) {
				force[i] += temp[i];
			}
		}
	}
	return force;
}

std::vector<double> Simulation::force(Particle a, Particle b) {
	std::vector<double> displacement = Simulation::displacement(a,b);
	double distance = 0;
	for (double component : Simulation::displacement(a, b)) {
		distance += std::pow(component, 2);
	}
	distance = std::sqrt(distance);
	double f = -24 * (2*std::pow(1/distance, 13) - std::pow(1/distance, 7));
	std::vector<double> force = std::vector<double>();
	for (double component : displacement) {
		force.push_back(component * f / distance);
	}
	return force;
} 

void Simulation::wrapParticle(Particle a) {
	a.setPosition(Simulation::positionInBox(a));
}