#include "sim.h"
#include <fstream>
#include <iostream>
#include "utils.h"
#include "potential.h"

using namespace MSM_MD_NS;
using namespace std;

double Sim::DELTA = 0.002;
double Sim::TEMP = 0.831;
vector<double> Sim::L = { 6.8, 6.8, 6.8 };
vector<bool> Sim::periodic = { true, true, true };

int Sim::timestep = 0;
int Sim::run_for = 0;
bool Sim::dumped = false;
bool Sim::thermoed = false;

double Sim::inv_t_damp_sq = 400;
double Sim::inv_t_set = 1.204;
double Sim::zeta = 0.0;

Atoms Sim::atoms = Atoms();
vector<double>(*Sim::force)(int) = Potential::lennard_jones_f_cutoff;
double (*Sim::potential)(int) = Potential::lennard_jones_e_cutoff;

vector<double>(*Sim::ij_force)(int, int) = Potential::lennard_jones_f_cutoff;
double (*Sim::ij_potential)(int, int) = Potential::lennard_jones_e_cutoff;

double (*Sim::scalar_force)(double) = Potential::lennard_jones_f_cutoff_scalar;

void (*Sim::integrator)() = Sim::verlet_one;

int Sim::dump_freq = 1;
int Sim::thermo_freq = 1;
int Sim::log_freq = 1;

string Sim::dumpfile = "";
string Sim::thermofile = "";


void Sim::run_sim(int n_steps)
{
	if (timestep == 0)
	{
		Logger::log("Time | Total E | Kinetic E | Potential E | Press | Virial Pressure | Energy Pressure | Temp \n---------------------------------------------");
	}

	for (int _ = 0; _ < n_steps; _++)
	{
		if (timestep % dump_freq == 0)
		{
			dump();
		}
		if (timestep % log_freq == 0)
		{
			log_out();
		}
		if (timestep % thermo_freq == 0)
		{
			thermo();
		}

		integrator();

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
		vector<double> fi = force(i);
		fx[i] = fi[0];
		fy[i] = fi[1];
		fz[i] = fi[2];
		//Logger::log(to_string(fx[i]));
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


		// PBCs
		if (periodic[0]) {
			if (atoms.x[i] < 0)
			{
				atoms.x[i] += L[0];
				atoms.x_flag[i] -= 1;
			}
			else if (atoms.x[i] > L[0])
			{
				atoms.x[i] -= L[0];
				atoms.x_flag[i] += 1;
			}
		}
		if (periodic[1]) {
			if (atoms.y[i] < 0)
			{
				atoms.y[i] += L[1];
				atoms.y_flag[i] -= 1;
			}
			else if (atoms.y[i] > L[1])
			{
				atoms.y[i] -= L[1];
				atoms.y_flag[i] += 1;
			}
		}
		if (periodic[2]) {
			if (atoms.z[i] < 0)
			{
				atoms.z[i] += L[2];
				atoms.z_flag[i] -= 1;
			}
			else if (atoms.z[i] > L[2])
			{
				atoms.z[i] -= L[2];
				atoms.z_flag[i] += 1;
			}
		}
	}

	// Calculate forces on new positions
	for (int i = 0; i < n; i++)
	{
		vector<double> fi = force(i);
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

void Sim::nose_hoover_one()
{
	int n = atoms.n_atoms;

	vector<double> fx(n);
	vector<double> fy(n);
	vector<double> fz(n);

	// Calculate forces on initial positions
	for (int i = 0; i < n; i++)
	{
		vector<double> fi = force(i);
		fx[i] = fi[0];
		fy[i] = fi[1];
		fz[i] = fi[2];
		//Logger::log(to_string(fx[i]));
	}

	// Calculate half-timestep velocities
	for (int i = 0; i < n; i++)
	{
		int t = atoms.type[i];

		atoms.vx[i] += (DELTA / 2) * (fx[i] / atoms.mass[t] - zeta * atoms.vx[i]);
		atoms.vy[i] += (DELTA / 2) * (fy[i] / atoms.mass[t] - zeta * atoms.vy[i]);
		atoms.vz[i] += (DELTA / 2) * (fz[i] / atoms.mass[t] - zeta * atoms.vz[i]);
	}

	// Calculate new positions
	for (int i = 0; i < n; i++)
	{
		atoms.x[i] += DELTA * atoms.vx[i];
		atoms.y[i] += DELTA * atoms.vy[i];
		atoms.z[i] += DELTA * atoms.vz[i];


		// PBCs
		if (periodic[0]) {
			if (atoms.x[i] < 0)
			{
				atoms.x[i] += L[0];
				atoms.x_flag[i] -= 1;
			}
			else if (atoms.x[i] > L[0])
			{
				atoms.x[i] -= L[0];
				atoms.x_flag[i] += 1;
			}
		}
		if (periodic[1]) {
			if (atoms.y[i] < 0)
			{
				atoms.y[i] += L[1];
				atoms.y_flag[i] -= 1;
			}
			else if (atoms.y[i] > L[1])
			{
				atoms.y[i] -= L[1];
				atoms.y_flag[i] += 1;
			}
		}
		if (periodic[2]) {
			if (atoms.z[i] < 0)
			{
				atoms.z[i] += L[2];
				atoms.z_flag[i] -= 1;
			}
			else if (atoms.z[i] > L[2])
			{
				atoms.z[i] -= L[2];
				atoms.z_flag[i] += 1;
			}
		}
	}

	// Calculate new zeta
	double temp = calc_t();
	zeta += DELTA * inv_t_damp_sq * (temp * inv_t_set - 1);


	// Calculate forces on new positions
	for (int i = 0; i < n; i++)
	{
		vector<double> fi = force(i);
		fx[i] = fi[0];
		fy[i] = fi[1];
		fz[i] = fi[2];
	}

	// Calculate new velocites
	double prefac = 1 / (1 + DELTA * 0.5 * zeta);
	for (int i = 0; i < n; i++)
	{
		int t = atoms.type[i];

		atoms.vx[i] = (atoms.vx[i] + (DELTA / 2) * fx[i] / atoms.mass[t]) * prefac;
		atoms.vy[i] = (atoms.vy[i] + (DELTA / 2) * fy[i] / atoms.mass[t]) * prefac;
		atoms.vz[i] = (atoms.vz[i] + (DELTA / 2) * fz[i] / atoms.mass[t]) * prefac;
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
		pe += potential(i);
		//Logger::log(to_string(potential(i)));
	}

	return pe / 2; // Divide by two because of double counting ji and ij bonds
}

double Sim::calc_e()
{
	return calc_ke() + calc_pe();
}

double Sim::calc_t()
{
	double ke = calc_ke();

	return 2 * ke / utils::kB / 3 / atoms.n_atoms;
}

double Sim::calc_press()
{
	return calc_press_ide() + calc_press_vir();
}

double Sim::calc_press_ide()
{
	double T = calc_t();

	double V = (L[0] * L[1] * L[2]);

	return atoms.n_atoms * utils::kB * T / V;
}

double Sim::calc_press_vir()
{
	double V = (L[0] * L[1] * L[2]);

	vector<double> x = atoms.x;
	vector<double> y = atoms.y;
	vector<double> z = atoms.z;

	double virial = 0.0;

	for (int i = 0; i < atoms.n_atoms; i++)
	{
		for (int j = i + 1; j < atoms.n_atoms; j++)
		{
			if (i == j)
			{
				continue;
			}

			double dx = periodic[0] ? utils::periodic_dist(x[j] - x[i], L[0]) : x[j] - x[i];
			double dy = periodic[1] ? utils::periodic_dist(y[j] - y[i], L[1]) : y[j] - y[i];
			double dz = periodic[2] ? utils::periodic_dist(z[j] - z[i], L[2]) : z[j] - z[i];

			double r = sqrt(dx * dx + dy * dy + dz * dz);
			double f = scalar_force(r);
			virial += r * f;
		}

		//vector<double> f = force(i);
		//virial += f[0] * atoms.x[i] + f[1] * atoms.y[i] + f[2] * atoms.z[i];
	}

	return virial / V / 3;
}

vector<double> Sim::calc_momentum()
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

double Sim::calc_msd()
{
	double msd = 0.0;

	for (int i = 0; i < atoms.n_atoms; i++)
	{
		double dx = atoms.x_flag[i] * L[0] + atoms.x[i] - atoms.x0[i];
		double dy = atoms.y_flag[i] * L[1] + atoms.y[i] - atoms.y0[i];
		double dz = atoms.z_flag[i] * L[2] + atoms.z[i] - atoms.z0[i];

		msd += (dx * dx + dy * dy + dz * dz);
	}

	return msd / atoms.n_atoms;
}

void Sim::dump()
{
	ofstream dump_io(dumpfile, ofstream::out | ofstream::app);
	if (dumped) {
		ofstream dump_io(dumpfile, ofstream::out | ofstream::app);
	}
	else {
		ofstream dump_io(dumpfile, ofstream::out | ofstream::trunc);
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

	dump_io << "ITEM: BOX BOUNDS xx yy zz" << endl;
	dump_io << "0.0 " << to_string(L[0]) << endl;
	dump_io << "0.0 " << to_string(L[1]) << endl;
	dump_io << "0.0 " << to_string(L[2]) << endl;

	// Write atom info
	dump_io << "ITEM: ATOMS id type x y z vx vy vz" << endl; // this whole method needs to be remade to be more generalized lol
	for (int i = 0; i < atoms.n_atoms; i++)
	{
		int id = atoms.id[i];
		int type = atoms.type[i];
		double x = atoms.x[i];
		double y = atoms.y[i];
		double z = atoms.z[i];
		double vx = atoms.vx[i];
		double vy = atoms.vy[i];
		double vz = atoms.vz[i];

		dump_io << id << " " << type << " " << to_string(x) << " " << to_string(y) << " " << to_string(z) << " " << to_string(vx) << " " << to_string(vy) << " " << to_string(vz) << endl;
	}

	dump_io.close();
}

void Sim::log_out()
{
	Logger::log(to_string(timestep * DELTA) + " " + to_string(Sim::calc_e()) + " " + to_string(Sim::calc_ke()) + " " + to_string(Sim::calc_pe()), true, true, "");

	double press = calc_press();
	double vir_press = calc_press_vir();
	double ke_press = calc_press_ide();
	double temp = calc_t();

	Logger::log(" " + to_string(press) + " " + to_string(vir_press) + " " + to_string(ke_press) + " " + to_string(temp));
}

void Sim::thermo()
{
	ofstream thermo_io(thermofile, ofstream::out | ofstream::app);
	if (thermoed) {
		ofstream thermo_io(thermofile, ofstream::out | ofstream::app);
	}
	else {
		ofstream thermo_io(thermofile, ofstream::out | ofstream::trunc);
		thermoed = true;

		thermo_io << "timestep time temperature pressure v_pressure e_pressure kinetic potential total mom_x mom_y mom_z msd" << endl;
	}

	if (!thermo_io.is_open())
	{
		return;
	}

	vector<double> p = calc_momentum();

	thermo_io << timestep << " " << timestep * DELTA << " " << calc_t() << " " << calc_press() << " " << calc_press_vir() << " " << calc_press_ide() << " " << calc_ke() << " " << calc_pe() << " " << calc_e() << " " << p[0] << " " << p[1] << " " << p[2] << " " << calc_msd() << endl;
	
	thermo_io.close();
}

void Sim::change_dump(string filename, int freq)
{
	dumped = false;
	dumpfile = filename;
	dump_freq = freq;
}

void Sim::change_thermo(string filename, int freq)
{
	thermoed = false;
	thermofile = filename;
	thermo_freq = freq;
}