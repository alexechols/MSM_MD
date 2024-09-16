#include "sim.h"
#include <fstream>
#include <iostream>
#include "utils.h"
#include "potential.h"

using namespace MSM_MD_NS;
using namespace std;

double Sim::DELTA = 0.002;
int Sim::timestep = 0;
int Sim::run_for = 0;
bool Sim::dumped = false;

Atoms Sim::atoms = Atoms();
vector<double>(*Sim::force)(Atoms,int) = Potential::lennard_jones_f;
double (*Sim::potential)(Atoms, int) = Potential::lennard_jones_e;

const char* Sim::dumpfile = "./msm_md.dump";


void Sim::verlet(int n_steps)
{
	Logger::log("Total E | Kinetic E | Potential E | Px | Py | Pz\n-----------------------------");
	for (int _ = 0; _ < n_steps; _++)
	{
		if (timestep % dump_freq == 0)
		{
			dump(dumpfile);
		}

		Logger::log(to_string(Sim::calc_e()) + " " + to_string(Sim::calc_ke()) + " " + to_string(Sim::calc_pe()), true, true, "");

		vector<double> p = Sim::calc_p();

		Logger::log(" " + to_string(p[0]) + " " + to_string(p[1]) + " " + to_string(p[2]));

		verlet_one();

		timestep++;
	}
}

void Sim::verlet_one()
{
	int n = atoms.n_atoms;

	vector<double> fx(n);
	vector<double> fy(n);
	vector<double> fz(n);

	// Calculate forces on initial positions
	for (int i = 0; i < n; i++)
	{
		vector<double> fi = force(atoms, i);
		fx[i] = fi[0];
		fy[i] = fi[1];
		fz[i] = fi[2];
	}

	// Calculate half-timestep velocities
	for (int i = 0; i < n; i++)
	{
		int t = atoms.type[i];

		atoms.vx[i] += (DELTA / 2) * fx[i] / atoms.mass[t];
		atoms.vy[i] += (DELTA / 2) * fy[i] / atoms.mass[t];
		atoms.vz[i] += (DELTA / 2) * fz[i] / atoms.mass[t];
	}

	// Calculate new positions
	for (int i = 0; i < n; i++)
	{
		atoms.x[i] += DELTA * atoms.vx[i];
		atoms.y[i] += DELTA * atoms.vy[i];
		atoms.z[i] += DELTA * atoms.vz[i];
	}

	// Calculate forces on new positions
	for (int i = 0; i < n; i++)
	{
		vector<double> fi = force(atoms, i);
		fx[i] = fi[0];
		fy[i] = fi[1];
		fz[i] = fi[2];
	}

	// Calculate new velocites
	for (int i = 0; i < n; i++)
	{
		int t = atoms.type[i];

		atoms.vx[i] += (DELTA / 2) * fx[i] / atoms.mass[t];
		atoms.vy[i] += (DELTA / 2) * fy[i] / atoms.mass[t];
		atoms.vz[i] += (DELTA / 2) * fz[i] / atoms.mass[t];
	}

}

double Sim::calc_ke()
{
	double ke = 0.0;

	vector<double> vx = atoms.vx;
	vector<double> vy = atoms.vy;
	vector<double> vz = atoms.vz;

	for (int i = 0; i < atoms.n_atoms; i++)
	{
		double v_mag_sq = vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];

		int t = atoms.type[i];

		ke += atoms.mass[t] * v_mag_sq;
	}

	return ke / 2; // Divide by two at the end for compute speed
}

double Sim::calc_pe()
{
	double pe = 0.0;

	for (int i = 0; i < atoms.n_atoms; i++)
	{
		pe += potential(atoms, i);
	}

	return pe / 2; // Divide by two because of double counting ji and ij bonds
}

double Sim::calc_e()
{
	return calc_ke() + calc_pe();
}

vector<double> Sim::calc_p()
{
	vector<double> p = { 0.0, 0.0, 0.0 };

	for (int i = 0; i < atoms.n_atoms; i++)
	{
		int t = atoms.type[i];
		double m = atoms.mass[t];

		p[0] += m * atoms.vx[i];
		p[1] += m * atoms.vy[i];
		p[2] += m * atoms.vz[i];
	}

	return p;
}

void Sim::dump(const char* filename)
{
	ofstream dump_io(filename, ofstream::out | ofstream::app);
	if (dumped) {
		ofstream dump_io(filename, ofstream::out | ofstream::app);
	}
	else {
		ofstream dump_io(filename, ofstream::out | ofstream::trunc);
		dumped = true;
	}
	
	if (!dump_io.is_open())
	{
		return;
	}

	// Roughly follows LAMMPS dump format (for ovito purposes)

	dump_io << "ITEM: TIMESTEP" << endl;
	dump_io << timestep << endl;

	dump_io << "ITEM: NUMBER OF ATOMS" << endl;
	dump_io << atoms.n_atoms << endl;

	// Write atom info
	dump_io << "ITEM: ATOMS id type x y z" << endl; // this whole method needs to be remade to be more generalized lol
	for (int i = 0; i < atoms.n_atoms; i++)
	{
		int id = atoms.id[i];
		int type = atoms.type[i];
		double x = atoms.x[i];
		double y = atoms.y[i];
		double z = atoms.z[i];

		dump_io << id << " " << type << " " << to_string(x) << " " << to_string(y) << " " << to_string(z) << endl;
	}

	dump_io.close();
}