#include "timer.h"
#include "logger.h"
#include <chrono>
#include <ctime>

using namespace std;
using namespace MSM_MD_NS;

Timer::Timer() {
	count = 0;
	elapsed = 0.0;
	running = false;
}

void Timer::start()
{
	if (!running)
	{
		start_time = chrono::system_clock::now();
		running = true;
	}
	else {
		Logger::warning("Cannot start timer while it is running!");
	}
}

void Timer::stop()
{
	if (running)
	{
		auto end = chrono::system_clock::now();

		chrono::duration<double> delta = end - start_time;

		elapsed += delta.count();
		count += 1;
		running = false;
	}
	else {
		Logger::warning("Cannot stop timer while it is not running!");
	}
}

double Timer::avg()
{
	if (count > 0) {
		return elapsed / count;
	}
	return 0;
}