#include "input.h";
#include <string>;
#include "sim.h"
#include "potential.h"
#include "logger.h"
#include "utils.h"

using namespace std;
using namespace MSM_MD_NS;

// Very nonexistent error checking rn

void Input::parse(int argc, char* argv[])
{
	string log_path = "./msm_md.log";
	int timesteps = 0;

	for (int i = 0; i < argc; i++)
	{
		string arg = argv[i];

		if (arg.compare("-in") == 0)
		{
			Sim::atoms = Atoms::create_atoms(argv[i + 1]); //Do this later so it gets properly logged (right now the logger stream is not initialized)
			Sim::force = Potential::lennard_jones_f_cutoff;
			Sim::potential = Potential::lennard_jones_e_cutoff;
		}

		else if (arg.compare("-log") == 0)
		{
			log_path = argv[i + 1];
		}
		else if (arg.compare("-dump") == 0)
		{
			Sim::dumpfile = argv[i + 1];
		}
		else if (arg.compare("-n") == 0)
		{
			Sim::run_for = utils::toInteger(argv[i + 1]);
		}
	}

	Logger::logfile = fopen(log_path.c_str(), "w");
}