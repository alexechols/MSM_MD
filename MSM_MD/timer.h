#ifndef TIMER
#define TIMER

#include <chrono>
#include <ctime>


using namespace std;

namespace MSM_MD_NS {
	class Timer {
	public:
		double elapsed;
		int count;

		Timer();

		void start();
		void stop();
		double avg();

	private:
		chrono::system_clock::time_point start_time;
		bool running;
	};
}

#endif