#include "potential.h"
#include <cmath>
#include "sim.h"
#include "utils.h"

using namespace std;
using namespace MSM_MD_NS;

double Potential::cutoff = 3.4;
bool Potential::cached = false;
double Potential::cutoff_e = 0.0;
double Potential::cutoff_e_deriv = 0.0;

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

		double dx = utils::periodic_dist(x[i] - x[j], Sim::L[0]);
		double dy = utils::periodic_dist(y[i] - y[j], Sim::L[1]);
		double dz = utils::periodic_dist(z[i] - z[j], Sim::L[2]);

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

		double dx = utils::periodic_dist(x[i] - x[j], Sim::L[0]);
		double dy = utils::periodic_dist(y[i] - y[j], Sim::L[1]);
		double dz = utils::periodic_dist(z[i] - z[j], Sim::L[2]);

		double inv_mag_sq = 1 / (dx * dx + dy * dy + dz * dz);
		double inv_mag_six = inv_mag_sq * inv_mag_sq * inv_mag_sq;

		e += 4 * (inv_mag_six * inv_mag_six - inv_mag_six);

	}

	return e;
}

vector<double> Potential::lennard_jones_f_cutoff(Atoms atoms, int i)
{
	if (!cached)
	{
		cache_cutoff();
	}

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

		double dx = utils::periodic_dist(x[i] - x[j], Sim::L[0]);
		double dy = utils::periodic_dist(y[i] - y[j], Sim::L[1]);
		double dz = utils::periodic_dist(z[i] - z[j], Sim::L[2]);

		//Logger::log(to_string(dx) + " " + to_string(dy) + " " + to_string(dz));

		double r_sq = dx * dx + dy * dy + dz * dz;
		double r = sqrt(r_sq);


		if (r <= cutoff)
		{
			double inv_mag_sq = 1 / r_sq;
			double inv_mag = 1 / r;
			double inv_mag_six = inv_mag_sq * inv_mag_sq * inv_mag_sq;

			double f_mag = 24 * (2 * inv_mag_six * inv_mag_six * inv_mag - inv_mag_six * inv_mag) + cutoff_e_deriv;

			//Logger::log(to_string(f_mag) + " " + to_string(r));

			force[0] += f_mag * dx * inv_mag;
			force[1] += f_mag * dy * inv_mag;
			force[2] += f_mag * dz * inv_mag;
		}

	}

	return force;

}

double Potential::lennard_jones_e_cutoff(Atoms atoms, int i) 
{
	if (!cached)
	{
		cache_cutoff();
	}

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

		double dx = utils::periodic_dist(x[i] - x[j], Sim::L[0]);
		double dy = utils::periodic_dist(y[i] - y[j], Sim::L[1]);
		double dz = utils::periodic_dist(z[i] - z[j], Sim::L[2]);

		double r_sq = dx * dx + dy * dy + dz * dz;
		double r = sqrt(r_sq);
		
		if (r <= cutoff)
		{
			double inv_mag_sq = 1 / r_sq;
			double inv_mag_six = inv_mag_sq * inv_mag_sq * inv_mag_sq;

			e += 4 * (inv_mag_six * inv_mag_six - inv_mag_six);
			e -= cutoff_e;
			e -= cutoff_e_deriv * (r - cutoff);
		}
		
	}

	return e;
}

void Potential::cache_cutoff()
{
	// Only LJ right now, need to generalize for arbitrary potentials later

	double inv_r = 1 / cutoff;
	double inv_six = inv_r * inv_r * inv_r * inv_r * inv_r * inv_r;

	cutoff_e = 4 * (inv_six * inv_six - inv_six);
	cutoff_e_deriv = 24 * (2 * inv_six * inv_six * inv_r - inv_six * inv_r);

	cached = true;
}