#include <vector>

#include "particle.h"
#include <fstream>

class Simulation {
	public:
		Simulation(int d, double l);
		std::vector<Particle> getParticles();
		std::vector<double> displacement(Particle, Particle);
		void initCubic(unsigned n, double temperature, double mass);
		double distance(Particle, Particle);
		std::vector<double> positionInBox(Particle&);
		double totalEnergy();
		double kineticEnergy();
		double potentialEnergy();
		double getLength();
		void updateParticles(double);
		std::vector<std::vector<double>> internalForces();
	private:
		int dim;
		double boxLength;
		std::vector<Particle> particles;
		double potential(Particle, Particle);
		std::vector<double> interaction(Particle, Particle);
		void wrapParticle(Particle&);
};