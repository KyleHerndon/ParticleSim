#include "particle.h"

#include <vector>

Particle::Particle(int dim, double m) {
	position = std::vector<double>(dim, 0);
	velocity = std::vector<double>(dim, 0);
	acceleration = std::vector<double>(dim, 0);
	mass = m;
}

std::vector<double> Particle::getPosition() {
	return position;
}

void Particle::setPosition(std::vector<double> pos) {
	position = pos;
}

std::vector<double> Particle::getVelocity() {
	return velocity;
}

void Particle::setVelocity(std::vector<double> vel) {
	velocity = vel;
}

std::vector<double> Particle::getAcceleration() {
	return acceleration;
}

void Particle::setAcceleration(std::vector<double> accel) {
	acceleration = accel;
}

double Particle::getMass() {
	return mass;
}