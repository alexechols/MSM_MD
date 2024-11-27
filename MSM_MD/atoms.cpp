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
}

void Atoms::add_atom(double x_pos, double y_pos, double z_pos, double vel_x, double vel_y, double vel_z, double m, double r)
{
	x.push_back(x_pos);
	y.push_back(y_pos);
	z.push_back(z_pos);

	x0.push_back(x_pos);
	y0.push_back(y_pos);
	z0.push_back(z_pos);

	fx.push_back(0.0);
	fy.push_back(0.0);
	fz.push_back(0.0);

	x_flag.push_back(0);
	y_flag.push_back(0);
	z_flag.push_back(0);

	vx.push_back(vel_x);
	vy.push_back(vel_y);
	vz.push_back(vel_z);

	px.push_back(vx[vx.size() - 1] * m);
	py.push_back(vy[vy.size() - 1] * m);
	pz.push_back(vz[vz.size() - 1] * m);

	id.push_back(id.size());
	mass.push_back(m);
	radius.push_back(r);

	n_atoms = id.size();
}

void Atoms::zero_momentum()
{
	double mu_x = 0.0;
	double mu_y = 0.0;
	double mu_z = 0.0;

	for (int i = 0; i < n_atoms; i++)
	{
		mu_x += vx[i] * mass[i];
		mu_y += vy[i] * mass[i];
		mu_z += vz[i] * mass[i];
	}

	mu_x /= n_atoms;
	mu_y /= n_atoms;
	mu_z /= n_atoms;

	double ke = 0.0;

	for (int i = 0; i < n_atoms; i++)
	{
		vx[i] -= mu_x / mass[i];
		vy[i] -= mu_y / mass[i];
		vz[i] -= mu_z / mass[i];

		px[i] = vx[i] * mass[i];
		py[i] = vy[i] * mass[i];
		pz[i] = vz[i] * mass[i];

		ke += mass[i] * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
	}
	Sim::ke = ke * 0.5;

}

void Atoms::reset_zero_positions()
{
	for (int i = 0; i < n_atoms; i++)
	{
		x0[i] = x[i];
		y0[i] = y[i];
		z0[i] = z[i];

		x_flag[i] = 0;
		y_flag[i] = 0;
		z_flag[i] = 0;
	}
}

void Atoms::read_add_atom(vector<string> header, vector<string> line)
{
	double x = 0.0, y = 0.0, z = 0.0;
	double vx = 0.0, vy = 0.0, vz = 0.0;
	double m = 1.0;
	double r = 1.0;
	//int id = -1;

	for (int i = 0; i < line.size(); i++)
	{
		if (i > header.size())
		{
			Logger::error("Too many parameters to parse in atomic data.");
		}

		string name = header[i];
		string el = line[i];

		if (!utils::isNumeric(el))
		{
			return;
		}

		if (name == "x") {
			x = utils::toFloat(el);
		}
		else if (name == "y")
		{
			y = utils::toFloat(el);
		}
		else if (name == "z")
		{
			z = utils::toFloat(el);
		}
		else if (name == "vx")
		{
			vx = utils::toFloat(el);
		}
		else if (name == "vy")
		{
			vy = utils::toFloat(el);
		}
		else if (name == "vz")
		{
			vz = utils::toFloat(el);
		}
		else if (name == "m")
		{
			m = utils::toFloat(el);
		}
		else if (name == "r")
		{
			r = utils::toFloat(el);
		}
		//else if (name == "id")
		//{
		//	if (!utils::isInteger(el))
		//	{
		//		Logger::error("Given id is not an integer.");
		//	}
		//	id = utils::toInteger(el);
		//}
		else {
			Logger::error("Unrecognized parameter \"" + name + "\" in header");
		}
	}

	add_atom(x, y, z, vx, vy, vz, m, r);
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
	int l_ind = 0;
	vector<string> header = {"x", "y", "z", "m", "vx", "vy", "vz", "r"};


	while (getline(data_io, line))
	{
		vector<string> split_line = utils::split(line);

		if (l_ind == 0)
		{
			// Check for header, if none is found, the file is assumed to contain data as follows
			// x y z mass vx vy vz r
			bool is_header = true;
			for (int i = 0; i < split_line.size(); i++)
			{
				if (utils::isNumeric(split_line[i]))
				{
					is_header = false;
					break;
				}
			}

			if (is_header)
			{
				header = split_line;
			}
			else {
				atoms.read_add_atom(header, split_line);
			}
		}

		else {
			atoms.read_add_atom(header, split_line);
		}

	}

	data_io.close();

	atoms.zero_momentum();

	Logger::log("Succesfully loaded ", true, true, "");
	Logger::log(to_string(atoms.n_atoms), true, true, "");
	Logger::log(" atoms\n");

	return atoms;
}

