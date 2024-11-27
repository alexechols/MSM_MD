#ifndef COMPUTE
#define COMPUTE

#include <vector>

using namespace std;
namespace MSM_MD_NS
{
	static class Computes
	{
	public: 
		static double calc_ke(); // Calculate kinetic energy at current moment
		static double calc_pe(); // Calculate potential energy at current moment
		static double calc_e(); // Calculate total energy at current moment
		static double calc_t(); // Instantaneous temperature
		static double calc_press(); // Instantaneous pressure
		static double calc_press_ide(); //Non-virial component of pressure
		static double calc_press_vir(); //Virial component of pressure
		static double calc_msd(); // Mean squared displacement
		static vector<double> calc_momentum(); // Calculates the net momentum in x, y, z
	};
}

#endif // !COMPUTE
