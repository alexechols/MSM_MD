#ifndef SIM
#define SIM

#include "atoms.h"
#include <vector>

using namespace std;
namespace MSM_MD_NS
{
	static class Sim {
	public:
		static double DELTA;
		static double TEMP;
		static vector<double> L;

		static int timestep;
		static int run_for;
		static const char* dumpfile;

		static const int dump_freq = 100;

		static Atoms atoms;
		static vector<double> (*force)(Atoms, int);
		static double (*potential)(Atoms, int);

		static void verlet(int n_steps); // Performs verlet integration over n_steps steps with length DELTA each

		static double calc_ke(); // Calculate kinetic energy at current moment
		static double calc_pe(); // Calculate potential energy at current moment
		static double calc_e(); // Calculate total energy at current moment
		static double calc_t(); // Instantaneous temperature
		static double calc_press(); // Instantaneous pressure
		static vector<double> calc_p(); // Calculates the net momentum in x, y, z

		static void dump(const char* filename);
		static void log_out();

	private:
		static void verlet_one();
		static bool dumped;
	};
}



#endif // !SIM
