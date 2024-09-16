#include "potential.h"
#include <cmath>

using namespace std;
using namespace MSM_MD_NS;

vector<double> Potential::lennard_jones_f(Atoms atoms, int i)
{
	vector<double> force = { 0.0, 0.0, 0.0 };

	vector<double> x = atoms.x;
	vector<double> y = atoms.y;
	vector<double> z = atoms.z;

	for (int j = 0; j < atoms.n_atoms; j++)
	{
		if (i == j)
		{
			continue; // Obviously don't compute forces with self
		}

		double dx = x[i] - x[j];
		double dy = y[i] - y[j];
		double dz = z[i] - z[j];

		double inv_mag_sq = 1 / (dx * dx + dy * dy + dz * dz);
		double inv_mag = sqrt(inv_mag_sq);
		double inv_mag_six = inv_mag_sq * inv_mag_sq * inv_mag_sq;

		double f_mag = 24 * (2 * inv_mag_six * inv_mag_six * inv_mag - inv_mag_six * inv_mag);

		force[0] += f_mag * dx * inv_mag;
		force[1] += f_mag * dy * inv_mag;
		force[2] += f_mag * dz * inv_mag;
	}

	return force;
}

double Potential::lennard_jones_e(Atoms atoms, int i)
{
	double e = 0;

	vector<double> x = atoms.x;
	vector<double> y = atoms.y;
	vector<double> z = atoms.z;

	for (int j = 0; j < atoms.n_atoms; j++)
	{
		if (i == j)
		{
			continue; // Obviously don't compute potential with self
		}

		double dx = x[i] - x[j];
		double dy = y[i] - y[j];
		double dz = z[i] - z[j];

		double inv_mag_sq = 1 / (dx * dx + dy * dy + dz * dz);
		double inv_mag_six = inv_mag_sq * inv_mag_sq * inv_mag_sq;

		e += 4 * (inv_mag_six * inv_mag_six - inv_mag_six);

	}

	return e;
}