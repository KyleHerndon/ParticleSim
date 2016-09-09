#include <vector>

#include "particle.h"

class Simulation {
	public:
		Simulation(double l);
		std::vector<double> getDisplacement(Particle a, Particle b);
		double getDistance(Particle a, Particle b);
		std::vector<double> positionInBox(Particle a);
	private:
		int dim;
		double boxLength;
		std::vector<Particle> particles;

}