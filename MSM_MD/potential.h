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

		static vector<double> lennard_jones_f(int i); // Calculates the force on the i-th atom using LJ
		static double lennard_jones_e(int i);

		static vector<double> lennard_jones_f(int i, int j); // Calculates the force on the i-th atom from the j-th atom using LJ
		static double lennard_jones_e(int i, int j);

		static double lennard_jones_f_scalar(double r);

		static vector<double> lennard_jones_f_cutoff(int i);
		static double lennard_jones_e_cutoff(int i);

		static vector<double> lennard_jones_f_cutoff(int i, int j);
		static double lennard_jones_e_cutoff(int i, int j);
		
		static double lennard_jones_f_cutoff_scalar(double r);

	private:
		static double cutoff_e; // Actual potential value at cutoff
		static double cutoff_e_deriv; // Actual e derivative value at cutoff
		static bool cached; // If cutoff_e and cutoff_f have been properly calculated

		static void cache_cutoff(); // Caches the cutoff values for e and f
	};
}

#endif