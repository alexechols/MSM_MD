#include "atoms.h"

#include <iostream>
#include <fstream>
#include "logger.h"
#include "utils.h"
#include "sim.h"
#include "random.h"
#include <math.h>

using namespace MSM_MD_NS;
using namespace std;

Atoms::Atoms()
{
	n_atoms = 0;
	mass.resize(1);
	mass[0] = 1.0;
}

void Atoms::add_atom(double px, double py, double pz)
{
	x.push_back(px);
	y.push_back(py);
	z.push_back(pz);

	vx.push_back(Random::gaussian(0.0, sqrt(Sim::TEMP)));
	vy.push_back(Random::gaussian(0.0, sqrt(Sim::TEMP)));
	vz.push_back(Random::gaussian(0.0, sqrt(Sim::TEMP)));

	type.push_back(0);
	id.push_back(id.size());

	n_atoms = id.size();
}

void Atoms::zero_momentum()
{
	double mu_x = 0.0;
	double mu_y = 0.0;
	double mu_z = 0.0;

	for (int i = 0; i < n_atoms; i++)
	{
		mu_x += vx[i];
		mu_y += vy[i];
		mu_z += vz[i];
	}

	mu_x /= n_atoms;
	mu_y /= n_atoms;
	mu_z /= n_atoms;

	for (int i = 0; i < n_atoms; i++)
	{
		vx[i] -= mu_x;
		vy[i] -= mu_y;
		vz[i] -= mu_z;
	}
}

Atoms Atoms::create_atoms(const char *filename)
{
	Logger::log("Reading atomic data from ", true, true, "");
	Logger::log(filename);

	ifstream data_io(filename);
	Atoms atoms = Atoms();

	if (!data_io.is_open())
	{
		Logger::error("Failed to load data");
	}

	string line;

	while (getline(data_io, line))
	{
		vector<string> split_line = utils::split(line);

		if (split_line.size() == 3)
		{
			// All lines of length 3 and all numeric data will be taken as atomic position data
			vector<double> p = { 0.0, 0.0, 0.0 };
			bool is_atomic = true;

			for (int i = 0; i < 3; i++)
			{
				if (!utils::isNumeric(split_line[i]))
				{
					is_atomic = false;
				}
				else {
					p[i] = utils::toFloat(split_line[i]);
				}
			}

			if (is_atomic)
			{
				atoms.add_atom(p[0], p[1], p[2]);
			}

		}
	}

	data_io.close();

	atoms.zero_momentum();

	Logger::log("Succesfully loaded ", true, true, "");
	Logger::log(to_string(atoms.n_atoms), true, true, "");
	Logger::log(" atoms\n");

	return atoms;
}