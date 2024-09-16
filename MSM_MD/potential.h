#ifndef POTENTIAL
#define POTENTIAL

#include <vector>
#include "atoms.h"

using namespace std;

namespace MSM_MD_NS
{
	static class Potential {
	public:
		static vector<double> lennard_jones_f(Atoms atoms, int i); // Calculates the force on the i-th atom using LJ
		static double lennard_jones_e(Atoms atoms, int i);
	};
}

#endif