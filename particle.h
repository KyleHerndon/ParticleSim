#include <vector>

class Particle {
	public:
		Particle();
		Particle(double x, double y, double z);
		std::vector<double> getPosition();
		std::vector<double> getVelocity();
		double getDistance(Particle other);
		double getMass();
	private:
		std::vector<double> position;
		std::vector<double> velocity;
		double mass;
};