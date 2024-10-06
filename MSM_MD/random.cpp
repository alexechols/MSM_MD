#include "random.h"
#include <numbers>
#include "logger.h"

#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;
using namespace MSM_MD_NS;

unsigned int Random::x = 123456789;
unsigned int Random::y = 362436069;
unsigned int Random::z = 521288629;
unsigned int Random::w = 88675123;

unsigned int Random::next() {
	int t = x ^ (x << 11);

	x = y;
	y = z;
	z = w;

	w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));

	return w;
}

double Random::uniform(double min, double max)
{
	double r = (double)(next());
	
	return (min - max) * uniform() + max;
}

double Random::uniform()
{
	return (float)(next() % 0x800000) / 0x7FFFFF;
}

double Random::gaussian(double mu, double sigma)
{
	double r1 = uniform(0, 1);
	double r2 = uniform(0, 1);

	double n = sqrt(-2 * log(r1)) * cos(2 * M_PI * r2);

	return mu + sigma * n;
}

void Random::seed(unsigned int seed)
{
	x = seed;
	y = (unsigned int)MAGIC * x + 1;
	z = (unsigned int)MAGIC * y + 1;
	w = (unsigned int)MAGIC * z + 1;
}