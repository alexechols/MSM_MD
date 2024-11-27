#include "sim.h"
#include <fstream>
#include <iostream>
#include "utils.h"
#include "potential.h"
#include "computes.h"
#include "constants.h"


using namespace MSM_MD_NS;
using namespace std;

double Sim::DELTA = 0.002;
vector<double> Sim::L = { 6.8, 6.8, 6.8 };
vector<bool> Sim::periodic = { false, false, false };

double Sim::virial = 0.0;
double Sim::pe = 0.0;
double Sim::ke = 0.0;
double Sim::mom = 0.0;

int Sim::timestep = 0;
bool Sim::dumped = false;
bool Sim::thermoed = false;

double Sim::inv_t_damp_sq = 400;
double Sim::inv_t_set = 1.204;
double Sim::zeta = 0.0;

Atoms Sim::atoms = Atoms();

void (*Sim::update_forces)() = Potential::lj_update_forces_potentials;
void (*Sim::integrator)() = Sim::verlet_one;

int Sim::dump_freq = 1;
int Sim::thermo_freq = 1;
int Sim::log_freq = 1;

string Sim::dumpfile = "";
string Sim::thermofile = "";

Timer Sim::force_timer = Timer();
Timer Sim::integrator_timer = Timer();
Timer Sim::io_timer = Timer();
Timer Sim::global = Timer();


void Sim::run_sim(int n_steps)
{
	if (timestep == 0)
	{
		Logger::log("Time | Total E | Kinetic E | Potential E | Press | Virial Pressure | Energy Pressure | Temp \n---------------------------------------------");
		update_forces_wrapper();
		dump();
		log_out();
		thermo();
	}

	for (int _ = 0; _ < n_steps; _++)
	{

		integrator_wrapper();

		timestep++;

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
	}
}

void Sim::update_forces_wrapper()
{
	force_timer.start();
	update_forces();
	force_timer.stop();
}

void Sim::integrator_wrapper()
{
	integrator_timer.start();
	integrator();
	integrator_timer.stop();
}

void Sim::verlet_one()
{
	int n = atoms.n_atoms;
	
	//update_forces();

	// Calculate half-timestep velocities
	for (int i = 0; i < n; i++)
	{

		atoms.vx[i] += (DELTA / 2) * atoms.fx[i] / atoms.mass[i];
		atoms.vy[i] += (DELTA / 2) * atoms.fy[i] / atoms.mass[i];
		atoms.vz[i] += (DELTA / 2) * atoms.fz[i] / atoms.mass[i];
		
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

	update_forces_wrapper();

	ke = 0.0;
	mom = 0.0;

	// Calculate new velocites
	for (int i = 0; i < n; i++)
	{	
		atoms.vx[i] += (DELTA / 2) * atoms.fx[i] / atoms.mass[i];
		atoms.vy[i] += (DELTA / 2) * atoms.fy[i] / atoms.mass[i];
		atoms.vz[i] += (DELTA / 2) * atoms.fz[i] / atoms.mass[i];

		atoms.px[i] = atoms.vx[i] * atoms.mass[i];
		atoms.py[i] = atoms.vy[i] * atoms.mass[i];
		atoms.pz[i] = atoms.vz[i] * atoms.mass[i];

		double ke_add = atoms.mass[i] * (atoms.vx[i] * atoms.vx[i] + atoms.vy[i] * atoms.vy[i] + atoms.vz[i] * atoms.vz[i]);
		ke += ke_add;
	}

	ke *= 0.5;
}

void Sim::nose_hoover_one()
{
	int n = atoms.n_atoms;

	//update_forces();

	// Calculate half-timestep velocities
	for (int i = 0; i < n; i++)
	{
		atoms.vx[i] += (DELTA / 2) * (atoms.fx[i] / atoms.mass[i] - zeta * atoms.vx[i]);
		atoms.vy[i] += (DELTA / 2) * (atoms.fy[i] / atoms.mass[i] - zeta * atoms.vy[i]);
		atoms.vz[i] += (DELTA / 2) * (atoms.fz[i] / atoms.mass[i] - zeta * atoms.vz[i]);
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
	double temp = Computes::calc_t();
	zeta += DELTA * inv_t_damp_sq * (temp * inv_t_set - 1);

	update_forces_wrapper();

	ke = 0.0;
	mom = 0.0;

	// Calculate new velocites
	double prefac = 1 / (1 + DELTA * 0.5 * zeta);
	for (int i = 0; i < n; i++)
	{

		atoms.vx[i] = (atoms.vx[i] + (DELTA / 2) * atoms.fx[i] / atoms.mass[i]) * prefac;
		atoms.vy[i] = (atoms.vy[i] + (DELTA / 2) * atoms.fy[i] / atoms.mass[i]) * prefac;
		atoms.vz[i] = (atoms.vz[i] + (DELTA / 2) * atoms.fz[i] / atoms.mass[i]) * prefac;

		atoms.px[i] = atoms.vx[i] * atoms.mass[i];
		atoms.py[i] = atoms.vy[i] * atoms.mass[i];
		atoms.pz[i] = atoms.vz[i] * atoms.mass[i];

		double ke_add = atoms.mass[i] * (atoms.vx[i] * atoms.vx[i] + atoms.vy[i] * atoms.vy[i] + atoms.vz[i] * atoms.vz[i]);
		ke += ke_add;
		mom += sqrt(ke_add * atoms.mass[i]);
	}

	ke *= 0.5;
}

void Sim::yoshida_one()
{
	int n = atoms.n_atoms;

	// Calculate x1 positions
	for (int i = 0; i < n; i++)
	{
		atoms.x[i] += DELTA * atoms.vx[i] * Const::c1;
		atoms.y[i] += DELTA * atoms.vy[i] * Const::c1;
		atoms.z[i] += DELTA * atoms.vz[i] * Const::c1;


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

	// Calculate v1 velocities
	update_forces_wrapper();

	for (int i = 0; i < n; i++)
	{
		atoms.vx[i] += Const::d1 * DELTA * atoms.fx[i] / atoms.mass[i];
		atoms.vy[i] += Const::d1 * DELTA * atoms.fy[i] / atoms.mass[i];
		atoms.vz[i] += Const::d1 * DELTA * atoms.fz[i] / atoms.mass[i];
	}

	// Calculate x2 positions
	for (int i = 0; i < n; i++)
	{
		atoms.x[i] += DELTA * atoms.vx[i] * Const::c2;
		atoms.y[i] += DELTA * atoms.vy[i] * Const::c2;
		atoms.z[i] += DELTA * atoms.vz[i] * Const::c2;


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


	// Calculate v2 velocities
	update_forces_wrapper();

	for (int i = 0; i < n; i++)
	{
		atoms.vx[i] += Const::d2 * DELTA * atoms.fx[i] / atoms.mass[i];
		atoms.vy[i] += Const::d2 * DELTA * atoms.fy[i] / atoms.mass[i];
		atoms.vz[i] += Const::d2 * DELTA * atoms.fz[i] / atoms.mass[i];
	}

	// Calculate x3 positions
	for (int i = 0; i < n; i++)
	{
		atoms.x[i] += DELTA * atoms.vx[i] * Const::c2;
		atoms.y[i] += DELTA * atoms.vy[i] * Const::c2;
		atoms.z[i] += DELTA * atoms.vz[i] * Const::c2;


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


	// Calculate v3 velocities (Final)
	update_forces_wrapper();
	ke = 0.0;
	mom = 0.0;

	for (int i = 0; i < n; i++)
	{
		atoms.vx[i] += Const::d1 * DELTA * atoms.fx[i] / atoms.mass[i];
		atoms.vy[i] += Const::d1 * DELTA * atoms.fy[i] / atoms.mass[i];
		atoms.vz[i] += Const::d1 * DELTA * atoms.fz[i] / atoms.mass[i];

		atoms.px[i] = atoms.vx[i] * atoms.mass[i];
		atoms.py[i] = atoms.vy[i] * atoms.mass[i];
		atoms.pz[i] = atoms.vz[i] * atoms.mass[i];

		double ke_add = atoms.mass[i] * (atoms.vx[i] * atoms.vx[i] + atoms.vy[i] * atoms.vy[i] + atoms.vz[i] * atoms.vz[i]);
		ke += ke_add;
	}

	ke *= 0.5;


	// Calculate x4 positions (Final)
	for (int i = 0; i < n; i++)
	{
		atoms.x[i] += DELTA * atoms.vx[i] * Const::c1;
		atoms.y[i] += DELTA * atoms.vy[i] * Const::c1;
		atoms.z[i] += DELTA * atoms.vz[i] * Const::c1;


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

}

void Sim::dump()
{
	io_timer.start();
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
	dump_io << "ITEM: ATOMS id x y z vx vy vz mass radius" << endl; // this whole method needs to be remade to be more generalized lol
	for (int i = 0; i < atoms.n_atoms; i++)
	{
		int id = atoms.id[i];
		double x = atoms.x[i];
		double y = atoms.y[i];
		double z = atoms.z[i];
		double vx = atoms.vx[i];
		double vy = atoms.vy[i];
		double vz = atoms.vz[i];
		double m = atoms.mass[i];
		double r = atoms.radius[i];

		dump_io << id << " " << to_string(x) << " " << to_string(y) << " " << to_string(z) << " " << to_string(vx) << " " << to_string(vy) << " " << to_string(vz) << " " << to_string(m) << " " << to_string(r) <<  endl;
	}

	dump_io.close();
	io_timer.stop();
}

void Sim::log_out()
{
	io_timer.start();
	Logger::log(to_string(timestep * DELTA) + " " + to_string(Computes::calc_e()) + " " + to_string(Computes::calc_ke()) + " " + to_string(Computes::calc_pe()), true, true, "");

	double press = Computes::calc_press();
	double vir_press = Computes::calc_press_vir();
	double ke_press = Computes::calc_press_ide();
	double temp = Computes::calc_t();

	Logger::log(" " + to_string(press) + " " + to_string(vir_press) + " " + to_string(ke_press) + " " + to_string(temp));
	io_timer.stop();
}

void Sim::thermo()
{
	io_timer.start();
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

	vector<double> p = Computes::calc_momentum();

	thermo_io << timestep << " " << timestep * DELTA << " " << Computes::calc_t() << " " << Computes::calc_press() << " " << Computes::calc_press_vir() << " " << Computes::calc_press_ide() << " " << Computes::calc_ke() << " " << Computes::calc_pe() << " " << Computes::calc_e() << " " << p[0] << " " << p[1] << " " << p[2] << " " << Computes::calc_msd() << endl;
	
	thermo_io.close();
	io_timer.stop();
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

void Sim::report_times()
{
	double tot = global.elapsed;
	double force = force_timer.elapsed;
	double io = io_timer.elapsed;
	double integrate = integrator_timer.elapsed - force;
	double other = tot - force - io - integrate;

	double p_force = force / tot * 100;
	double p_io = io / tot * 100;
	double p_integrate = integrate / tot * 100;
	double p_other = 100 - p_force - p_io - p_integrate;

	Logger::log("Elapsed Time: " + to_string(tot) + "s");
	double dt = tot / Sim::timestep;
	Logger::log("Avg Time per Timestep: " + to_string(dt) + "s");
	Logger::log("\nBreakdown:\n\tForce Calculations: " + to_string(force) + "s (" + to_string(p_force) + "%)");
	Logger::log("\tIntegrator (Excluding Forces): " + to_string(integrate) + "s (" + to_string(p_integrate) + "%)");
	Logger::log("\tI/O (Dump/Log/Thermo): " + to_string(io) + "s (" + to_string(p_io) + "%)");
	Logger::log("\tOther: " + to_string(other) + "s (" + to_string(p_other) + "%)");
}