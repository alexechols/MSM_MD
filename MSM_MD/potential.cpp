#include "potential.h"
#include <cmath>
#include "sim.h"
#include "utils.h"
#include "logger.h"

using namespace std;
using namespace MSM_MD_NS;

double Potential::cutoff = 2.5;
double Potential::cutoff_sq = 2.5 * 2.5;

bool Potential::cached = false;

double Potential::cutoff_e = 0.0;
double Potential::cutoff_e_deriv = 0.0;

void Potential::lj_update_forces_potentials()
{
	Atoms atoms = Sim::atoms;

	int n = atoms.n_atoms;
	vector<double> x = atoms.x;
	vector<double> y = atoms.y;
	vector<double> z = atoms.z;

	for (int i = 0; i < n; i++)
	{
		atoms.fx[i] = 0.0;
		atoms.fy[i] = 0.0;
		atoms.fz[i] = 0.0;
	}


	double virial = 0.0;
	double pe = 0.0;

	for (int i = 0; i < n; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			double dx = Sim::periodic[0] ? utils::periodic_dist(x[j] - x[i], Sim::L[0]) : x[j] - x[i];
			double dy = Sim::periodic[1] ? utils::periodic_dist(y[j] - y[i], Sim::L[1]) : y[j] - y[i];
			double dz = Sim::periodic[2] ? utils::periodic_dist(z[j] - z[i], Sim::L[2]) : z[j] - z[i];

			double r_two = dx * dx + dy * dy + dz * dz;
			double r = sqrt(r_two);
			double inv_r_two = 1 / r_two;
			double inv_r_six = inv_r_two * inv_r_two * inv_r_two;

			double inv_r = 1 / r;

			double force_factor = -24 * (2 * inv_r_six - 1) * inv_r_six * inv_r;
			double pot_factor = 4 * (inv_r_six - 1) * inv_r_six;

			pe += pot_factor;
			virial -= force_factor * r;

			atoms.fx[i] += force_factor * dx * inv_r;
			atoms.fy[i] += force_factor * dy * inv_r;
			atoms.fz[i] += force_factor * dz * inv_r;

			atoms.fx[j] -= force_factor * dx * inv_r;
			atoms.fy[j] -= force_factor * dy * inv_r;
			atoms.fz[j] -= force_factor * dz * inv_r;
		}
	}

	Sim::virial = virial;
	Sim::pe = pe;
	Sim::atoms = atoms;
}

void Potential::lj_cut_update_forces_potentials()
{
	if (!cached)
	{
		cache_cutoff();
	}

	//AtoSim::atoms;

	int n = Sim::atoms.n_atoms;
	vector<double> x = Sim::atoms.x;
	vector<double> y = Sim::atoms.y;
	vector<double> z = Sim::atoms.z;

	for (int i = 0; i < n; i++)
	{
		Sim::atoms.fx[i] = 0;
		Sim::atoms.fy[i] = 0;
		Sim::atoms.fz[i] = 0;
	}

	double virial = 0.0;
	double pe = 0.0;

	for (int i = 0; i < n; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			double dx = Sim::periodic[0] ? utils::periodic_dist(x[j] - x[i], Sim::L[0]) : x[j] - x[i];
			double dy = Sim::periodic[1] ? utils::periodic_dist(y[j] - y[i], Sim::L[1]) : y[j] - y[i];
			double dz = Sim::periodic[2] ? utils::periodic_dist(z[j] - z[i], Sim::L[2]) : z[j] - z[i];

			double r_two = dx * dx + dy * dy + dz * dz;

			if (r_two > cutoff_sq)
			{
				continue;
			}

			double r = sqrt(r_two);
			double inv_r_two = 1 / r_two;
			double inv_r_six = inv_r_two * inv_r_two * inv_r_two;

			double inv_r = 1 / r;

			double force_factor = -24 * (2 * inv_r_six - 1) * inv_r_six * inv_r + cutoff_e_deriv;
			double pot_factor = 4 * (inv_r_six - 1) * inv_r_six - cutoff_e + (r - cutoff) * cutoff_e_deriv;

			pe += pot_factor;
			virial -= force_factor * r;

			//Logger::log(to_string(r) + " " + to_string(force_factor));

			Sim::atoms.fx[i] += force_factor * dx * inv_r;
			Sim::atoms.fy[i] += force_factor * dy * inv_r;
			Sim::atoms.fz[i] += force_factor * dz * inv_r;

			Sim::atoms.fx[j] -= force_factor * dx * inv_r;
			Sim::atoms.fy[j] -= force_factor * dy * inv_r;
			Sim::atoms.fz[j] -= force_factor * dz * inv_r;
		}
	}


	Sim::virial = virial;
	Sim::pe = pe;
}

void Potential::gravity_update_forces_potentials()
{
	Atoms atoms = Sim::atoms;

	int n = atoms.n_atoms;
	vector<double> x = atoms.x;
	vector<double> y = atoms.y;
	vector<double> z = atoms.z;
	vector<double> m = atoms.mass;

	for (int i = 0; i < n; i++)
	{
		atoms.fx[i] = 0.0;
		atoms.fy[i] = 0.0;
		atoms.fz[i] = 0.0;
	}


	double virial = 0.0;
	double pe = 0.0;

	for (int i = 0; i < n; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			double dx = Sim::periodic[0] ? utils::periodic_dist(x[j] - x[i], Sim::L[0]) : x[j] - x[i];
			double dy = Sim::periodic[1] ? utils::periodic_dist(y[j] - y[i], Sim::L[1]) : y[j] - y[i];
			double dz = Sim::periodic[2] ? utils::periodic_dist(z[j] - z[i], Sim::L[2]) : z[j] - z[i];

			double r_two = dx * dx + dy * dy + dz * dz;
			double inv_r_two = 1 / r_two;
			double inv_r = sqrt(inv_r_two);
			double r = sqrt(r_two);

			double force_factor = m[i] * m[j] * inv_r_two;

			double pot_factor = -m[i] * m[j] * inv_r;

			pe += pot_factor;
			virial -= force_factor * r;

			atoms.fx[i] += force_factor * dx * inv_r;
			atoms.fy[i] += force_factor * dy * inv_r;
			atoms.fz[i] += force_factor * dz * inv_r;

			atoms.fx[j] -= force_factor * dx * inv_r;
			atoms.fy[j] -= force_factor * dy * inv_r;
			atoms.fz[j] -= force_factor * dz * inv_r;
		}
	}


	Sim::atoms = atoms;
	Sim::virial = virial;
	Sim::pe = pe;
}

void Potential::dummy_update_forces_potentials()
{
	Atoms atoms = Sim::atoms;

	int n = atoms.n_atoms;

	for (int i = 0; i < n; i++)
	{
		atoms.fx[i] = 0.0;
		atoms.fy[i] = 0.0;
		atoms.fz[i] = 0.0;
	}

	Sim::atoms = atoms;
	Sim::virial = 0.0;
	Sim::pe = 0.0;
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