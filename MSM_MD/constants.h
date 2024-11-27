#ifndef CONST
#define CONST

namespace MSM_MD_NS
{
	static class Const {
	public:
		static const double kB; // Boltzmann's constant
		static const double Na; // Avagadro's 

		static const double eps;
		static const double sigma;

		static const double w0; // Intermediate constant for yoshida integrator
		static const double w1; // Intermediate constant for yoshida integrator
		static const double c1; // Constant for yoshida integrator
		static const double c2; // Constant for yoshida integrator
		static const double d1; // Constant for yoshida integrator
		static const double d2; // Constant for yoshida integrator

	};
}

#endif