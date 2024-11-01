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

		vector<double> x0;
		vector<double> y0;
		vector<double> z0;

		vector<int> x_flag;
		vector<int> y_flag;
		vector<int> z_flag;

		vector<double> vx;
		vector<double> vy;
		vector<double> vz;

		vector<double> fx;
		vector<double> fy;
		vector<double> fz;

		vector<double> mass;

		int n_atoms;

		Atoms();

		void add_atom(double px, double py, double pz);

		static Atoms create_atoms(const char *filename);

		void zero_momentum();

		void reset_zero_positions();
		

	};
}

#endif // !ATOMS



