#ifndef SIM
#define SIM

#include "atoms.h"
#include <vector>
#include <string>
#include "timer.h"

using namespace std;
namespace MSM_MD_NS
{
	static class Sim {
	public:
		static double DELTA; // Timestep

		static vector<double> L; // Box dimensions
		static vector<bool> periodic; // PBCs

		static double virial; // Virial component of the pressure
		static double pe; // Potential Energy
		static double ke; // Kinetic Energy
		static double mom; // System Momementum

		static int timestep; // Current timestep
		static double time; // Current time (timestep * DELTA)

		// Parameters for nvt, not used otherwise
		static double inv_t_damp_sq;
		static double inv_t_set;
		static double zeta;

		// File paths
		static string dumpfile;
		static string thermofile;

		// Logging freq
		static int dump_freq;
		static int log_freq;
		static int thermo_freq;

		static Timer force_timer;
		static Timer integrator_timer;
		static Timer io_timer;
		static Timer global;

		static Atoms atoms;

		// Pointers for force update
		static void(*update_forces)();
		
		//Pointer for integrator
		static void (*integrator)();

		static void run_sim(int n_steps); // Performs verlet integration over n_steps steps with length DELTA each

		static void dump();
		static void thermo();
		static void log_out();
		static void report_times();

		static void change_dump(string filename, int freq);
		static void change_thermo(string filename, int freq);
		
		// Integrators
		static void verlet_one();
		static void yoshida_one();
		static void nose_hoover_one();

	private:
		static bool dumped;
		static bool thermoed;

		static void update_forces_wrapper(); // Wrapper for timing purposes
		static void integrator_wrapper(); // Wrapper for timing purposes
	};
}



#endif // !SIM
