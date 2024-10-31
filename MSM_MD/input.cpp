#include "input.h";
#include <string>;
#include "sim.h"
#include "potential.h"
#include "logger.h"
#include "utils.h"
#include <iostream>
#include <fstream>
#include "command.h"

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
			ifstream input_io(argv[i+1]);

			if (!input_io.is_open())
			{
				Logger::error("Failed to load input file");
			}

			string line;

			while (getline(input_io, line))
			{
				Command::parse_line(line);
			}

		}
		/*else
		{
			Logger::error("Invalid command line option \"" + (string) argv[i] + "\"");
		}*/
	}

	//Logger::logfile = fopen(log_path.c_str(), "w");
}