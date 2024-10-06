#ifndef RANDOM
#define RANDOM



namespace MSM_MD_NS {
	static class Random {
	public:
		static double gaussian(double mu, double sigma);
		static double uniform(double min, double max);
		static double uniform();
		static void seed(unsigned int seed);

	private:
		static unsigned int next();

		static unsigned int x;
		static unsigned int y;
		static unsigned int z;
		static unsigned int w;

		static const unsigned int MAGIC = 1812433253;
	};
}

#endif // !
