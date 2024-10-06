#ifndef ATOMS
#define ATOMS

#include <vector>

using namespace std;

namespace MSM_MD_NS {

	class Atoms {
	public:
		vector<int> id;
		vector<int> type;

		vector<double> x;
		vector<double> y;
		vector<double> z;

		vector<double> vx;
		vector<double> vy;
		vector<double> vz;

		vector<double> mass;

		int n_atoms;

		Atoms();

		void add_atom(double px, double py, double pz);

		static Atoms create_atoms(const char *filename);

	private:
		void zero_momentum();

	};
}

#endif // !ATOMS



