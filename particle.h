#include <valarray>

class Particle {
	public:
		Particle(int dim, double m);
		std::vector<double> getPosition();
		void setPosition(std::vector<double> pos);
		std::vector<double> getVelocity();
		void setVelocity(std::vector<double> vel);
		std::vector<double> getAcceleration();
		void setAcceleration(std::vector<double> accel);
		double getMass();
		double getEnergy();
	private:
		std::vector<double> position;
		std::vector<double> velocity;
		std::vector<double> acceleration;
		double mass;
};