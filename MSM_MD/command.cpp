#include "command.h"
#include "utils.h"
#include "sim.h"
#include "potential.h"
#include "atoms.h"
#include "random.h"
#include "constants.h"

using namespace std;
using namespace MSM_MD_NS;

void Command::parse_line(string line)
{
	vector<string> split_line = utils::split(line);

	if (split_line.size() < 1) {
		return;
	}

	if (split_line.size() == 1 && split_line[0] == "")
	{
		return;
	}

	string command = split_line[0];

	if (command == "log") { Command::log(line); }
	else if (command == "dump") { Command::dump(line); }
	else if (command == "thermo") { Command::thermo(line); }
	else if (command == "timestep") { Command::timestep(line); }
	else if (command == "potential") { Command::potential(line); }
	else if (command == "box") { Command::box(line); }
	else if (command == "periodic") { Command::periodic(line); }
	else if (command == "atoms") { Command::atoms(line); }
	else if (command == "velocity") { Command::velocity(line); }
	else if (command == "run") { Command::run(line); }
	else if (command == "seed") { Command::seed(line); }
	else {
		Logger::error("Malformed or unrecognized command \"" + split_line[0] + "\" in line \"" + line + "\"");
	}
}

void Command::log(string line)
{
	vector<string> split_line = utils::split(line);

	if (split_line.size() != 3)
	{
		Logger::error("Incorrect parameters in command \"" + line + "\"");
	}

	if (!utils::isInteger(split_line[2])) {
		Logger::error("Frequency should be an integer, value \"" + split_line[2] + "\" is invalid");
	}

	Logger::log("Switching log file to " + split_line[1] + "...");
	Logger::log("Switching log frequency to " + split_line[2] + " timesteps...");

	Logger::change_file(split_line[1]);
	Sim::log_freq = utils::toInteger(split_line[2]);
}

void Command::dump(string line)
{
	vector<string> split_line = utils::split(line);

	if (split_line.size() != 3)
	{
		Logger::error("Incorrect parameters in command \"" + line + "\"");
	}

	if (!utils::isInteger(split_line[2])) {
		Logger::error("Frequency should be an integer, value \"" + split_line[2] + "\" is invalid");
	}

	Logger::log("Switching dump file to " + split_line[1] + "...");
	Logger::log("Switching dump frequency to " + split_line[2] + " timesteps...");

	Sim::change_dump(split_line[1], utils::toInteger(split_line[2]));
}

void Command::thermo(string line)
{
	vector<string> split_line = utils::split(line);

	if (split_line.size() != 3)
	{
		Logger::error("Incorrect parameters in command \"" + line + "\"");
	}

	if (!utils::isInteger(split_line[2])) {
		Logger::error("Frequency should be an integer, value \"" + split_line[2] + "\" is invalid");
	}

	Logger::log("Switching thermo file to " + split_line[1] + "...");
	Logger::log("Switching thermo frequency to " + split_line[2] + " timesteps...");

	Sim::change_thermo(split_line[1], utils::toInteger(split_line[2]));
}

void Command::timestep(string line)
{
	vector<string> split_line = utils::split(line);

	if (split_line.size() != 2)
	{
		Logger::error("Incorrect parameters in command \"" + line + "\"");
	}

	if (!utils::isNumeric(split_line[1]))
	{
		Logger::error("Timestep should be a float, value \"" + split_line[1] + "\" is invalid");
	}

	Logger::log("Timestep: " + split_line[1]);
	Sim::DELTA = utils::toFloat(split_line[1].c_str());
}

void Command::potential(string line)
{
	vector<string> split_line = utils::split(line);

	if (split_line.size() != 2 && split_line.size() != 3)
	{
		Logger::error("Incorrect parameters in command \"" + line + "\"");
	}

	if (split_line.size() == 3 && !utils::isNumeric(split_line[2]))
	{
		Logger::error("Cutoff should be a float, value \"" + split_line[2] + "\" is invalid");
	}

	if (split_line[1] == "lennard_jones")
	{
		if (split_line.size() == 2)
		{
			Sim::force = Potential::lennard_jones_f;
			Sim::potential = Potential::lennard_jones_e;
			Sim::ij_force = Potential::lennard_jones_f;
			Sim::ij_potential = Potential::lennard_jones_e;
			Sim::scalar_force = Potential::lennard_jones_f_scalar;

			Logger::log("Potential: Lennard Jones");
		}
		else
		{
			Sim::force = Potential::lennard_jones_f_cutoff;
			Sim::potential = Potential::lennard_jones_e_cutoff;
			Sim::ij_force = Potential::lennard_jones_f_cutoff;
			Sim::ij_potential = Potential::lennard_jones_e_cutoff;
			Sim::scalar_force = Potential::lennard_jones_f_cutoff_scalar;

			double cut = utils::toFloat(split_line[2]);
			Potential::cutoff = cut;
			Potential::cutoff_sq = cut * cut;

			Logger::log("Potential: Lennard Jones (Cutoff: " + split_line[2] + ")");
		}
	}

	else {
		Logger::error("\"" + split_line[1] + "\" is not a valid potential");
	}
}

void Command::box(string line)
{
	vector<string> split_line = utils::split(line);

	if (split_line.size() != 4)
	{
		Logger::error("Incorrect parameters in command \"" + line + "\"");
	}

	vector<double> dims = { 0.0, 0.0, 0.0 };

	for (int i = 1; i < 4; i++)
	{
		if (!utils::isNumeric(split_line[i]))
		{
			Logger::error("Dimensions should be floats, value \"" + split_line[i] + "\" is invalid");
		}
		else
		{
			dims[i-1] = utils::toFloat(split_line[i]);
		}
	}

	Sim::L = dims;
	Logger::log("Setting box dimensions");
}

void Command::periodic(string line)
{
	vector<string> split_line = utils::split(line);

	if (split_line.size() != 2)
	{
		Logger::error("Incorrect parameters in command \"" + line + "\"");
	}

	vector<bool> periodic = { false, false, false };

	for (int i = 0; i < split_line[1].size(); i++)
	{
		char c =  split_line[1].at(i);

		if (c == 'x' || c == 'X')
		{
			periodic[0] = true;
		}
		else if (c == 'y' || c == 'Y')
		{
			periodic[1] = true;
		}
		else if (c == 'z' || c == 'Z')
		{
			periodic[2] = true;
		}
		else
		{
			string message = "Axis \"";
			message += c;
			Logger::error(message + "\" is not defined and cannot be made periodic");
		}
	}

	Sim::periodic = periodic;
	Logger::log("Setting boundary conditions");
}

void Command::atoms(string line)
{
	vector<string> split_line = utils::split(line);

	if (split_line.size() != 2)
	{
		Logger::error("Incorrect parameters in command \"" + line + "\"");
	}

	Sim::atoms = Atoms::create_atoms(split_line[1].c_str());
}

void Command::velocity(string line)
{
	vector<string> split_line = utils::split(line);

	if (split_line.size() != 2 && split_line.size() != 3)
	{
		Logger::error("Incorrect parameters in command \"" + line + "\"");
	}

	if (Sim::atoms.n_atoms == 0)
	{
		Logger::warning("Attempting to set velocity for zero atoms. Nothing will occur");
		return;
	}

	if (split_line.size() == 2)
	{
		if (split_line[1] == "zero")
		{
			for (int i = 0; i < Sim::atoms.n_atoms; i++)
			{
				Sim::atoms.vx[i] = 0.0;
				Sim::atoms.vy[i] = 0.0;
				Sim::atoms.vz[i] = 0.0;				
			}
			Logger::log("Velocities zeroed");
		}
		else {
			Logger::error("Unrecognized velocity mode \"" + split_line[1] + "\"");
		}
	}
	else if (split_line.size() == 3)
	{
		if (split_line[1] == "temp")
		{
			if (!utils::isNumeric(split_line[2]))
			{
				Logger::error("Temperature must be a float, value \"" + split_line[2] + "\" is invalid");
			}

			double temp = utils::toFloat(split_line[2]);

			for (int i = 0; i < Sim::atoms.n_atoms; i++)
			{
				Sim::atoms.vx[i] = Random::gaussian(0.0, sqrt(temp * Const::kB / Const::eps));
				Sim::atoms.vy[i] = Random::gaussian(0.0, sqrt(temp * Const::kB / Const::eps));
				Sim::atoms.vz[i] = Random::gaussian(0.0, sqrt(temp * Const::kB / Const::eps));
			}
			Sim::atoms.zero_momentum();
			Logger::log("Velocities initialized at T = " + split_line[2]);
		}
		else {
			Logger::error("Unrecognized velocity mode \"" + split_line[1] + "\"");
		}
	}
}

void Command::run(string line)
{
	vector<string> split_line = utils::split(line);

	if (split_line.size() != 2)
	{
		Logger::error("Incorrect parameters in command \"" + line + "\"");
	}

	if (!utils::isInteger(split_line[1]))
	{
		Logger::error("Timesteps should be an integer, value \"" + split_line[1] + "\" is invalid");
	}

	Logger::log("Running for " + split_line[1] + " timesteps");

	Sim::run_for = utils::toInteger(split_line[1]);
	Sim::run_sim(Sim::run_for);
}

void Command::seed(string line)
{
	vector<string> split_line = utils::split(line);

	if (split_line.size() != 2)
	{
		Logger::error("Incorrect parameters in command \"" + line + "\"");
	}

	if (!utils::isInteger(split_line[1]))
	{
		Logger::error("Seed should be an integer, value \"" + split_line[1] + "\" is invalid");
	}

	Random::seed(utils::toInteger(split_line[1]));
	Logger::log("Random generator seeded with " + split_line[1]);
}

void Command::integrate(string line)
{
	vector<string> split_line = utils::split(line);

	if (split_line.size() != 2 || split_line.size() != 4)
	{
		Logger::error("Incorrect parameters in command \"" + line + "\"");
	}

	string mode = split_line[1];

	if (mode == "nve")
	{
		if (split_line.size() != 2)
		{
			Logger::error("Mode nve does not take additional parameters");
		}

		Sim::integrator = Sim::verlet_one;
		Logger::log("Integrator is NVE");
	}
	else if (mode == "nvt")
	{
		if (split_line.size() != 4)
		{
			Logger::error("Mode nvt takes a temperature and a damping constant");
		}

		if (!utils::isNumeric(split_line[2]))
		{
			Logger::error("Temperature must be a float, value \"" + split_line[2] + "\" is invalid");
		}

		if (!utils::isNumeric(split_line[3]))
		{
			Logger::error("Damping constant must be a float, value \"" + split_line[3] + "\" is invalid");
		}

		Sim::t_damp = utils::toFloat(split_line[3]);
		Sim::t_set = utils::toFloat(split_line[2]);
		Sim::integrator = Sim::nose_hoover_one;
	}
}