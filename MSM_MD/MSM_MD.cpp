// MSM_MD.cpp : Defines the entry point for the application.
//

#include "MSM_MD.h"

using namespace std;
using namespace MSM_MD_NS;

int main(int argc, char* argv[])
{
	Input::parse(argc, argv);
	Sim::verlet(Sim::run_for);

	fclose(Logger::logfile);

}
