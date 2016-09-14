#include "simulation.h"

#include <cmath>
#include <math.h>
#include <iostream>
#include <fstream>
#include <ctime>


Simulation::Simulation(int d, double l) {
	dim = d;
	boxLength = l;
}

std::vector<Particle> Simulation::getParticles(){
	return particles;
}

void Simulation::initCubic(unsigned n, double temperature, double mass) {
	unsigned count = pow(n, 3);
	particles.reserve(count);
	double spacing = boxLength/n;
	double offset = (n % 2 == 0 ) ? boxLength/(2*n) : 0;
	offset -= spacing * floor(n/2);
	std::vector<double> centerOfMass (dim, 0);
	srand(time(NULL)); // randomize seed
	for (unsigned i = 0; i < n; i++) {
		for (unsigned j = 0; j < n; j++) {
			for (unsigned k = 0; k < n; k++) {
				Particle p (3, mass);
				p.setPosition(std::vector<double> 
					({i * spacing + offset, j * spacing + offset, k * spacing + offset}));
				p.setVelocity(std::vector<double>(
					{((double) rand())/RAND_MAX-.5, ((double) rand())/RAND_MAX-.5, ((double) rand())/RAND_MAX-.5}));
				for (int c = 0; c < dim; c++) {
					centerOfMass[c] += p.getVelocity()[c];

				}
				particles.push_back(p);
			}
		}
	}
	double total = 0;
	for (int c = 0; c < dim; c++) {
		centerOfMass[c] /= count;
	}
	for (unsigned i = 0; i < particles.size(); i++) {
		std::vector<double> v = particles[i].getVelocity();
		for (int c = 0; c < dim; c++) {
			v[c] -= centerOfMass[c];
			total += pow(v[c], 2);	
		}
		particles[i].setVelocity(v);
	}
	for (unsigned i = 0; i < particles.size(); i++) {
		std::vector<double> v = particles[i].getVelocity();
		for (int c = 0; c < dim; c++) {
			v[c] = v[c] * sqrt(3 * temperature * count / (mass * total));

		}
		particles[i].setVelocity(v);
	}
	// Set acceleration
	std::vector<std::vector<double>> forces = Simulation::internalForces();
	for (unsigned i = 0; i < particles.size(); i++) {
		std::vector<double> acc = forces[i];
		for (int c =0; c < dim; c++) {
			acc[c] /= mass;
		}
		particles[i].setAcceleration(acc);
	}
}

std::vector<double> Simulation::displacement(Particle a, Particle b) {
	std::vector<double> v = std::vector<double>(dim, 0);
	for (int i = 0; i < dim; i++) {
		double temp = b.getPosition()[i] - a.getPosition()[i];
		if (temp < (boxLength / 2) && temp > -(boxLength / 2)) {
			v[i] = temp;
		} else {
			v[i] = boxLength-temp;
			v[i] += boxLength/2;
			v[i] = fmod(v[i], boxLength);
			v[i] -= boxLength/2;
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

std::vector<double> Simulation::positionInBox(Particle &a) {
	std::vector<double> newPos = a.getPosition();
	for (unsigned i = 0; i < newPos.size(); i++) {
		newPos[i] += boxLength / 2;
		newPos[i] = (newPos[i] < 0) ? boxLength-newPos[i] : newPos[i];
		newPos[i] = fmod(newPos[i], boxLength);
		newPos[i] -= boxLength / 2;
	}
	return newPos;
}

double Simulation::totalEnergy() {
	return Simulation::kineticEnergy()+Simulation::potentialEnergy();
}

double Simulation::kineticEnergy() {
	double ret = 0;
	for (Particle &particle : particles) {
		ret += particle.getEnergy();
	}
	std::cout << "kinetic: " << ret << "\n";
	return ret;
}

double Simulation::potentialEnergy() {
	double ret = 0;
	for (unsigned i = 0; i < particles.size(); i++) {
		for (unsigned j = i +1; j < particles.size(); j++) {
			ret+=Simulation::potential(particles[i], particles[j]);
		}
	}
	std::cout <<  "potential: " << ret << "\n";
	return ret;
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
			newAcceleration[component] = forces[i][component] / particles[i].getMass();
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

double Simulation::potential(Particle a, Particle b) {
	std::vector<double> displacement = Simulation::displacement(a,b);
	double distance = 0;
	for (double component : displacement) {
		distance += std::pow(component, 2);
	}
	distance = std::sqrt(distance);
	
	return 4 * (std::pow(1/distance, 12) - std::pow(1/distance, 6));
}

std::vector<double> Simulation::interaction(Particle a, Particle b) {
	std::vector<double> displacement = Simulation::displacement(a, b);
	double distance = Simulation::distance(a, b);
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