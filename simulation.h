#include <vector>

#include "particle.h"

class Simulation {
	public:
		Simulation(int d, double l);
		std::vector<double> displacement(Particle a, Particle b);
		double distance(Particle a, Particle b);
		std::vector<double> positionInBox(Particle a);
		std::vector<double> internalForce(Particle a);
	private:
		int dim;
		double boxLength;
		std::vector<Particle> particles;
		std::vector<double> force(Particle a, Particle b);
		void wrapParticle(Particle a);
};