#ifndef COMMAND
#define COMMAND

#include <string>

using namespace std;


namespace MSM_MD_NS {
	class Command {
	public:
		static void parse_line(string line);

		static void log(string line);
		static void dump(string line);
		static void thermo(string line);
		static void timestep(string line);
		static void potential(string line);
		static void box(string line);
		static void periodic(string line);
		static void atoms(string line);
		static void velocity(string line);
		static void run(string line);
		static void seed(string line);
		static void integrate(string line);
		
		/* Commands 
		*	log [filepath] [int freq]: Change the path of the log file, frequency to log thermo values on
		*	dump [filepath] [int freq] : Change the path of the dump file, frequency to dump on
		*	thermo [filepath] [int freq] : Change the path of the thermo file, frequency to dump thermo on
		* 
		*	timestep [float timestep] : floating point value of timestep length. Default is 0.002
		*	potential [string potential] : string name of the potential to use. Passing a floating point value after the cutoff will implement a cutoff at that distance
		*	
		*	box [float x] [float y] [float z] : x,y,z dimensions of the simulation bounding box
		*	periodic [string mode] : string (ex. xy, xyz) for the dimensions to apply periodic boundary conditions
		*	atoms [filepath] : Load atomic data from filepath. Currently only loads positions and zeros velocities
		*	velocity [string mode] : Intialize atomic velocities according to some scheme
		*	
		*	run [int n] : run simulation for n timesteps
		*	integrate [string mode] : specifies integration mode (nve, nvt)
		* 
		*	seed [int x] : Seed the random generator with x
		*/
	};
}

#endif // !COMMAND