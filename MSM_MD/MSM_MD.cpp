// MSM_MD.cpp : Defines the entry point for the application.
//

#include "MSM_MD.h"

using namespace std;
using namespace MSM_MD_NS;

int main(int argc, char* argv[])
{
	Random::seed(1234567);
	

	Sim::global.start();

	Input::parse(argc, argv);

	Sim::global.stop();

	Sim::report_times();

	fclose(Logger::logfile);

}
