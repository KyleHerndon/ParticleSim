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

void Simulation::initCubic(unsigned n, double mass, double temperature) {
	unsigned count = pow(n, 3);
	particles.reserve(count);
	double spacing = boxLength/n;
	double offset = (n % 2 == 0 ) ? boxLength/(2*n) : 0;
	std::vector<double> centerOfMass (dim, 0);
	for (unsigned i = 0; i < n; i++) {
		for (unsigned j = 0; j < n; j++) {
			for (unsigned k = 0; k < n; k++) {
				Particle p (3, mass);
				p.setPosition(std::vector<double> 
					({i * spacing + offset, j * spacing + offset, k * spacing + offset}));
				p.setVelocity(std::vector<double>({rand()-.5, rand()-.5, rand()-.5}));
				for (int c = 0; c < dim; c++) {
					centerOfMass[c] += p.getVelocity()[c];
				}
				particles.push_back(p);
			}
		}
	}
	double total = 0;
	for (Particle particle : particles) {
		std::vector<double> v (dim, 0);
		for (int c = 0; c < dim; c++) {
			v[c] = particle.getVelocity()[c] - centerOfMass[c];
			total += pow(particle.getVelocity()[c], 2);
			particle.setVelocity(v);
		}
	}
	total = sqrt(total);
	for (Particle particle : particles) {
		std::vector<double> v (dim,0);
		for (int c = 0; c < dim; c++) {
			v[c] = 3 * particle.getVelocity()[c] * count * temperature / (mass * total) ;
		}
		particle.setVelocity(v);
	}
}

double Simulation::distance(Particle a, Particle b) {
	double ret = 0;
	for (double component : displacement(a, b)) {
		ret += std::pow(component, 2);
	}
	return std::sqrt(ret);
}

std::vector<double> Simulation::positionInBox(Particle &a) {
	std::vector<double> newPos = std::vector<double>(dim, 0);
	std::vector<double> pos = a.getPosition();
	for (unsigned i = 0; i < pos.size(); i++) {
		newPos[i] = (pos[i] > (boxLength / 2) ? boxLength-pos[i] : pos[i]);
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

double Simulation::getLength() {
	return boxLength;
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
		Simulation::wrapParticle(particles[i]);
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

void Simulation::addParticle(Particle a){

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