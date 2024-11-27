#ifndef POTENTIAL
#define POTENTIAL

#include <vector>
#include "atoms.h"

using namespace std;

namespace MSM_MD_NS
{
	static class Potential {
	public:
		static double cutoff;
		static double cutoff_sq;

		static void lj_update_forces_potentials();
		static void lj_cut_update_forces_potentials();
		static void gravity_update_forces_potentials();
		static void dummy_update_forces_potentials();

	private:
		static double cutoff_e; // Actual potential value at cutoff
		static double cutoff_e_deriv; // Actual e derivative value at cutoff
		static bool cached; // If cutoff_e and cutoff_f have been properly calculated

		static void cache_cutoff(); // Caches the cutoff values for e and f
	};
}

#endif