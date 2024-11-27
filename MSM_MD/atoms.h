#ifndef ATOMS
#define ATOMS

#include <vector>
#include <string>

using namespace std;

namespace MSM_MD_NS {

	class Atoms {
	public:
		vector<int> id;

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

		vector<double> px;
		vector<double> py;
		vector<double> pz;

		vector<double> fx;
		vector<double> fy;
		vector<double> fz;

		vector<double> mass;
		vector<double> radius;

		int n_atoms;

		Atoms();

		void add_atom(double x_pos, double y_pos, double z_pos, double vel_x, double vel_y, double vel_z, double m, double r);

		static Atoms create_atoms(const char *filename);
		void read_add_atom(vector<string> header, vector<string> line);

		void zero_momentum();

		void reset_zero_positions();
		

	};
}

#endif // !ATOMS



