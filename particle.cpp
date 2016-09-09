#include "particle.h"

#include <vector>

Particle::Particle(void) {
	position = std::vector<double>();
	velocity = std::vector<double>();
	mass = 0;
}

Particle::Particle(double x, double y, double z) {
	Particle();
	position.assign(0, x);
	position.assign(1, y);
	position.assign(2, z);
}

std::vector<double> Particle::getPosition() {
	return position;
}

std::vector<double> Particle::getVelocity() {
	return velocity;
}

double Particle::getMass() {
	return mass;
}

double Particle::getDistance(Particle other) {
	double ret = 0;
	for (unsigned i = 0; i < position.size(); i++) {
		ret += (position[i]*position[i]+other.position[i]*other.position[i]);
	}
	return ret;
}