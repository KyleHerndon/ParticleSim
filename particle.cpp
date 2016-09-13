#include "particle.h"

#include <vector>
#include <cmath>

Particle::Particle(int dim, double m) {
	position = std::vector<double>(dim, 0);
	velocity = std::vector<double>(dim, 0);
	acceleration = std::vector<double>(dim, 0);
	mass = m;
}

const std::vector<double> Particle::getPosition() {
	return position;
}

void Particle::setPosition(std::vector<double> pos) {
	position = pos;
}

const std::vector<double> Particle::getVelocity() {
	return velocity;
}

void Particle::setVelocity(std::vector<double> vel) {
	velocity = vel;
}

const std::vector<double> Particle::getAcceleration() {
	return acceleration;
}

void Particle::setAcceleration(std::vector<double> accel) {
	acceleration = accel;
}

double Particle::getMass() {
	return mass;
}

double Particle::getEnergy() {
	double energy = 0;
	for(double component : velocity) {
		energy += std::pow(component, 2);
	}
	energy *= mass * .5;
	return energy;
}