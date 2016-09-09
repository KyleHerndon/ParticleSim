#include <vector>

#include "particle.h"

class Simulation {
	public:
		Simulation(int d, double l);
		std::vector<double> displacement(Particle a, Particle b);
		double distance(Particle a, Particle b);
		std::vector<double> positionInBox(Particle a);
		double totalEnergy();
		double kineticEnergy();
		double potentialEnergy();
		void updateParticles(double timestep);
		std::vector<std::vector<double>> internalForces();
	private:
		int dim;
		double boxLength;
		std::vector<Particle> particles;
		double potential(Particle a, Particle b);
		std::vector<double> interaction(Particle a, Particle b);
		void wrapParticle(Particle a);
};