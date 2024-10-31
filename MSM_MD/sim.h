#ifndef SIM
#define SIM

#include "atoms.h"
#include <vector>
#include <string>

using namespace std;
namespace MSM_MD_NS
{
	static class Sim {
	public:
		static double DELTA;
		static double TEMP;
		static vector<double> L;
		static vector<bool> periodic;

		static int timestep;
		static int run_for;

		// Parameters for nvt, not used otherwise
		static double t_damp;
		static double t_set;

		static string dumpfile;
		static string thermofile;

		static int dump_freq;
		static int log_freq;
		static int thermo_freq;

		static Atoms atoms;

		// Pointers for force functions
		static vector<double> (*force)(int);
		static double (*potential)(int);
		static vector<double>(*ij_force)(int, int);
		static double(*ij_potential)(int, int);
		static double (*scalar_force)(double);
		
		//Pointer for integrator
		static void (*integrator)();

		static void run_sim(int n_steps); // Performs verlet integration over n_steps steps with length DELTA each

		static double calc_ke(); // Calculate kinetic energy at current moment
		static double calc_pe(); // Calculate potential energy at current moment
		static double calc_e(); // Calculate total energy at current moment
		static double calc_t(); // Instantaneous temperature
		static double calc_press(); // Instantaneous pressure
		static double calc_press_ide(); //Non-virial component of pressure
		static double calc_press_vir(); //Virial component of pressure
		static vector<double> calc_momentum(); // Calculates the net momentum in x, y, z

		static void dump();
		static void thermo();
		static void log_out();

		static void change_dump(string filename, int freq);
		static void change_thermo(string filename, int freq);
		
		// Integrators
		static void verlet_one();
		static void nose_hoover_one();
	private:
		static bool dumped;
		static bool thermoed;
	};
}



#endif // !SIM
