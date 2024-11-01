#include "potential.h"
#include <cmath>
#include "sim.h"
#include "utils.h"
#include "logger.h"

using namespace std;
using namespace MSM_MD_NS;

double Potential::cutoff = 3.4;
double Potential::cutoff_sq = 3.4;
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

	vector<double> fx = atoms.fx;
	vector<double> fy = atoms.fy;
	vector<double> fz = atoms.fz;

	for (int i = 0; i < n; i++)
	{
		fx[i] = 0;
		fy[i] = 0;
		fz[i] = 0;
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

			pe += 2 * pot_factor; // Times two because we are adding potential for ij and ji
			virial += force_factor * r;

			fx[i] += force_factor * dx * inv_r;
			fy[i] += force_factor * dy * inv_r;
			fz[i] += force_factor * dz * inv_r;

			fx[j] -= force_factor * dx * inv_r;
			fy[j] -= force_factor * dy * inv_r;
			fz[j] -= force_factor * dz * inv_r;
		}
	}

	Sim::virial = virial;
	Sim::pe = pe;
}

void Potential::lj_cut_update_force_potentials()
{
	if (!cached)
	{
		cache_cutoff();
	}

	Atoms atoms = Sim::atoms;

	int n = atoms.n_atoms;
	vector<double> x = atoms.x;
	vector<double> y = atoms.y;
	vector<double> z = atoms.z;

	vector<double> fx = atoms.fx;
	vector<double> fy = atoms.fy;
	vector<double> fz = atoms.fz;

	for (int i = 0; i < n; i++)
	{
		fx[i] = 0;
		fy[i] = 0;
		fz[i] = 0;
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

			double force_factor = -24 * (2 * inv_r_six - 1) * inv_r_six * inv_r - cutoff_e_deriv;
			double pot_factor = 4 * (inv_r_six - 1) * inv_r_six - cutoff_e + (r - cutoff) * cutoff_e_deriv;

			pe += 2 * pot_factor; // Times two because we are adding potential for ij and ji
			virial += force_factor * r;

			fx[i] += force_factor * dx * inv_r;
			fy[i] += force_factor * dy * inv_r;
			fz[i] += force_factor * dz * inv_r;

			fx[j] -= force_factor * dx * inv_r;
			fy[j] -= force_factor * dy * inv_r;
			fz[j] -= force_factor * dz * inv_r;
		}
	}

	Sim::virial = virial;
	Sim::pe = pe;
}

vector<double> Potential::lennard_jones_f(int i, int j)
{
	vector<double> force = { 0.0, 0.0, 0.0 };

	vector<double> x = Sim::atoms.x;
	vector<double> y = Sim::atoms.y;
	vector<double> z = Sim::atoms.z;

	if (i == j)
	{
		return force;
	}

	double dx = Sim::periodic[0] ? utils::periodic_dist(x[j] - x[i], Sim::L[0]) : x[j] - x[i];
	double dy = Sim::periodic[1] ? utils::periodic_dist(y[j] - y[i], Sim::L[1]) : y[j] - y[i];
	double dz = Sim::periodic[2] ? utils::periodic_dist(z[j] - z[i], Sim::L[2]) : z[j] - z[i];

	double inv_mag_sq = 1 / (dx * dx + dy * dy + dz * dz);
	//double inv_mag = sqrt(inv_mag_sq);
	double inv_mag_six = inv_mag_sq * inv_mag_sq * inv_mag_sq;

	//double f_mag = - 24 * (2 * inv_mag_six * inv_mag_six * inv_mag - inv_mag_six * inv_mag);
	double f_mag = -24 * (2 * inv_mag_six - 1) * inv_mag_six * inv_mag_sq;

	//force[0] = f_mag * dx * inv_mag;
	//force[1] = f_mag * dy * inv_mag;
	//force[2] = f_mag * dz * inv_mag;

	force[0] = f_mag * dx;
	force[1] = f_mag * dy;
	force[2] = f_mag * dz;

	//Logger::log(to_string(1 / inv_mag) + " " + to_string(f_mag));

	return force;
}

vector<double> Potential::lennard_jones_f(int i)
{
	vector<double> force = { 0.0, 0.0, 0.0 };

	for (int j = 0; j < Sim::atoms.n_atoms; j++)
	{
		vector<double> ij_force = lennard_jones_f(i, j);
		
		force[0] += ij_force[0];
		force[1] += ij_force[1];
		force[2] += ij_force[2];
	}

	return force;
}

double Potential::lennard_jones_f_scalar(double r)
{
	double inv_r = 1 / r;
	double inv_r_sq = inv_r * inv_r;
	double inv_r_six = inv_r_sq * inv_r_sq * inv_r_sq;

	return -24 * (2 * inv_r_six - 1) * inv_r_six * inv_r_sq;
}

double Potential::lennard_jones_e(int i, int j)
{
	double e = 0.0;

	if (i == j)
	{
		return e;
	}
	
	vector<double> x = Sim::atoms.x;
	vector<double> y = Sim::atoms.y;
	vector<double> z = Sim::atoms.z;

	double dx = Sim::periodic[0] ? utils::periodic_dist(x[j] - x[i], Sim::L[0]) : x[j] - x[i];
	double dy = Sim::periodic[1] ? utils::periodic_dist(y[j] - y[i], Sim::L[1]) : y[j] - y[i];
	double dz = Sim::periodic[2] ? utils::periodic_dist(z[j] - z[i], Sim::L[2]) : z[j] - z[i];

	double inv_mag_sq = 1 / (dx * dx + dy * dy + dz * dz);
	double inv_mag_six = inv_mag_sq * inv_mag_sq * inv_mag_sq;

	e = 4 * (inv_mag_six - 1) * inv_mag_six;
	
	return e;
}

double Potential::lennard_jones_e(int i)
{
	double e = 0;

	for (int j = 0; j < Sim::atoms.n_atoms; j++)
	{
		double ij_e = lennard_jones_e(i, j);

		e += ij_e;
	}

	return e;
}

vector<double> Potential::lennard_jones_f_cutoff(int i, int j)
{
	if (!cached)
	{
		cache_cutoff();
	}

	vector<double> force = { 0.0, 0.0, 0.0 };

	vector<double> x = Sim::atoms.x;
	vector<double> y = Sim::atoms.y;
	vector<double> z = Sim::atoms.z;

	if (i == j)
	{
		return force;
	}

	double dx = Sim::periodic[0] ? utils::periodic_dist(x[j] - x[i], Sim::L[0]) : x[j] - x[i];
	double dy = Sim::periodic[1] ? utils::periodic_dist(y[j] - y[i], Sim::L[1]) : y[j] - y[i];
	double dz = Sim::periodic[2] ? utils::periodic_dist(z[j] - z[i], Sim::L[2]) : z[j] - z[i];

	double r_sq = dx * dx + dy * dy + dz * dz;


	if (r_sq > cutoff_sq)
	{
		return force;
	}

	double inv_mag_sq = 1 / r_sq;
	double inv_mag = sqrt(inv_mag_sq);
	double inv_mag_six = inv_mag_sq * inv_mag_sq * inv_mag_sq;

	double f_mag = -24 * (2 * inv_mag_six * inv_mag_six * inv_mag - inv_mag_six * inv_mag) + cutoff_e_deriv;

	force[0] = f_mag * dx * inv_mag;
	force[1] = f_mag * dy * inv_mag;
	force[2] = f_mag * dz * inv_mag;

	return force;
}

vector<double> Potential::lennard_jones_f_cutoff(int i)
{

	vector<double> force = { 0.0, 0.0, 0.0 };

	for (int j = 0; j < Sim::atoms.n_atoms; j++)
	{
		vector<double> ij_force = lennard_jones_f_cutoff(i, j);

		force[0] += ij_force[0];
		force[1] += ij_force[1];
		force[2] += ij_force[2];
	}

	return force;

}

double Potential::lennard_jones_f_cutoff_scalar(double r)
{
	if (!cached)
	{
		cache_cutoff();
	}

	if (r > cutoff)
	{
		return 0.0;
	}

	double inv_r = 1 / r;
	double inv_r_sq = inv_r * inv_r;
	double inv_r_six = inv_r_sq * inv_r_sq * inv_r_sq;

	double f = 24 * (2 * inv_r_six - 1) * inv_r_six * inv_r;

	return f - cutoff_e_deriv;
}

double Potential::lennard_jones_e_cutoff(int i, int j)
{
	if (!cached)
	{
		cache_cutoff();
	}

	double e = 0;

	vector<double> x = Sim::atoms.x;
	vector<double> y = Sim::atoms.y;
	vector<double> z = Sim::atoms.z;

	if (i == j)
	{
		return e;
	}

	double dx = Sim::periodic[0] ? utils::periodic_dist(x[j] - x[i], Sim::L[0]) : x[j] - x[i];
	double dy = Sim::periodic[1] ? utils::periodic_dist(y[j] - y[i], Sim::L[1]) : y[j] - y[i];
	double dz = Sim::periodic[2] ? utils::periodic_dist(z[j] - z[i], Sim::L[2]) : z[j] - z[i];

	double r_sq = dx * dx + dy * dy + dz * dz;
	double r = sqrt(r_sq);

	if (r > cutoff)
	{
		return e;
	}

	double inv_mag_sq = 1 / r_sq;
	double inv_mag_six = inv_mag_sq * inv_mag_sq * inv_mag_sq;

	e += 4 * (inv_mag_six * inv_mag_six - inv_mag_six);
	e -= cutoff_e;
	e += cutoff_e_deriv * (r - cutoff);

	return e;
}
double Potential::lennard_jones_e_cutoff(int i) 
{
	double e = 0;

	for (int j = 0; j < Sim::atoms.n_atoms; j++)
	{
		e += lennard_jones_e_cutoff(i, j);
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