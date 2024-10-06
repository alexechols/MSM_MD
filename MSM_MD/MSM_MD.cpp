// MSM_MD.cpp : Defines the entry point for the application.
//

#include "MSM_MD.h"

using namespace std;
using namespace MSM_MD_NS;

int main(int argc, char* argv[])
{
	Random::seed(1234567);
	Input::parse(argc, argv);

	auto start = chrono::system_clock::now();

	Sim::verlet(Sim::run_for);

	auto end = chrono::system_clock::now();

	chrono::duration<double> elapsed = end - start;

	Logger::log("Elapsed Time: " + to_string(elapsed.count()) + "s");
	double dt = elapsed.count() / Sim::run_for;
	Logger::log("Avg Time per Timestep: " + to_string(dt) + "s");

	fclose(Logger::logfile);

}
