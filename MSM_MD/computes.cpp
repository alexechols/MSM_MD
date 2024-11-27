#include "computes.h"
#include "sim.h"
#include "utils.h"

using namespace std;
using namespace MSM_MD_NS;

double Computes::calc_ke()
{
	return Sim::ke;
}

double Computes::calc_pe()
{
	return Sim::pe;
}

double Computes::calc_e()
{
	return calc_ke() + calc_pe();
}

double Computes::calc_t()
{
	double ke = calc_ke();

	return 2 * ke / utils::kB / 3 / Sim::atoms.n_atoms;
}

double Computes::calc_press()
{
	return calc_press_ide() + calc_press_vir();
}

double Computes::calc_press_ide()
{
	vector<double> L = Sim::L;

	double T = calc_t();

	double V = (L[0] * L[1] * L[2]);

	return Sim::atoms.n_atoms * utils::kB * T / V;
}

double Computes::calc_press_vir()
{
	vector<double> L = Sim::L;

	double V = (L[0] * L[1] * L[2]);

	return Sim::virial / V / 3;
}

vector<double> Computes::calc_momentum()
{
	vector<double> p = { 0.0, 0.0, 0.0 };

	for (int i = 0; i < Sim::atoms.n_atoms; i++)
	{
		p[0] += Sim::atoms.px[i];
		p[1] += Sim::atoms.py[i];
		p[2] += Sim::atoms.pz[i];
	}

	return p;
}

double Computes::calc_msd()
{
	Atoms atoms = Sim::atoms;
	double msd = 0.0;

	for (int i = 0; i < atoms.n_atoms; i++)
	{
		double dx = atoms.x_flag[i] * Sim::L[0] + atoms.x[i] - atoms.x0[i];
		double dy = atoms.y_flag[i] * Sim::L[1] + atoms.y[i] - atoms.y0[i];
		double dz = atoms.z_flag[i] * Sim::L[2] + atoms.z[i] - atoms.z0[i];

		msd += (dx * dx + dy * dy + dz * dz);
	}

	return msd / atoms.n_atoms;
}
