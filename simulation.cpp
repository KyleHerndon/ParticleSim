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

std::vector<double> Simulation::positionInBox(const Particle &a) {
	std::vector<double> newPos = std::vector<double>();
	for (double component : a.getPosition()) {
		newPos.push_back((component > (boxLength / 2) ? boxLength-component : component));
	}
	return newPos;
}

double Simulation::totalEnergy() {
	return Simulation::kineticEnergy()+Simulation::potentialEnergy();
}

double Simulation::kineticEnergy() {
	double ret = 0;
	for (Particle particle : particles) {
		ret += particle.getEnergy();
	}
	return ret;
}

double Simulation::potentialEnergy() {
	double ret = 0;
	for (unsigned i = 0; i < particles.size(); i++) {
		for (unsigned j = i +1; j < particles.size(); j++) {
			ret+=Simulation::potential(particles[i], particles[j]);
		}
	}
	return 2*ret;
}

void Simulation::updateParticles(double timestep) {
	std::vector<std::vector<double>> forces = Simulation::internalForces();
	for (unsigned i = 0; i < particles.size(); i++) {
		std::vector<double> newPosition (dim, 0);
		std::vector<double> newVelocity (dim, 0);
		std::vector<double> newAcceleration (dim, 0);
		std::vector<double> oldPosition = particles[i].getPosition();
		std::vector<double> oldVelocity = particles[i].getVelocity();
		std::vector<double> oldAcceleration = particles[i].getAcceleration();
		for (int component = 0; component < dim; component++) {
			newAcceleration[component] /= particles[i].getMass();
			// Verlet Position
			newPosition[component] = oldPosition[component] + oldVelocity[component]*timestep+
				.5 * oldAcceleration[component] * pow(timestep, 2);
			// Verlet Velocity
			newVelocity[component] = oldVelocity[component] + 
				.5 * (oldAcceleration[component] + newAcceleration[component]) * timestep;
		}
		particles[i].setPosition(newPosition);
		particles[i].setVelocity(newVelocity);
		particles[i].setAcceleration(newAcceleration);
	}
}

std::vector<std::vector<double>> Simulation::internalForces() {
	// calculating all the particle forces at once, so I can get a 2x speed up using
	// Newton's 3rd Law
	std::vector<std::vector<double>> forces = std::vector<std::vector<double>>(
		particles.size(),std::vector<double>(dim, 0));
	for (unsigned i = 0; i < particles.size(); i++) {
		for (unsigned j = i+1; j < particles.size(); j++) {
			std::vector<double> temp = Simulation::interaction(particles[i], particles[j]);
			for (int component = 0; component < dim; component++) {
				forces[i][component] += temp[component];
				forces[j][component] -= temp[component];
			}
		}
	}
	return forces;
}

double Simulation::potential(Particle a, Particle b) {
	std::vector<double> displacement = Simulation::displacement(a,b);
	double distance = 0;
	for (double component : Simulation::displacement(a, b)) {
		distance += std::pow(component, 2);
	}
	distance = std::sqrt(distance);
	return 4 * (std::pow(1/distance, 12) + std::pow(1/distance, 6));
}

std::vector<double> Simulation::interaction(Particle a, Particle b) {
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

void Simulation::wrapParticle(Particle &a) {
	a.setPosition(Simulation::positionInBox(a));
}