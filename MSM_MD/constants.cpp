#include "constants.h"
#include "math.h"

using namespace MSM_MD_NS;

const double Const::kB = 1.380649E-23; // J / K
const double Const::Na = 6.0221415E23; // 1 / mol

const double Const::eps = 1.66E-21;
const double Const::sigma = 3.4;

const double Const::w0 = - cbrt(2) / (2 - cbrt(2));
const double Const::w1 = 1 / (2 - cbrt(2));
const double Const::c1 = w1 / 2; // C1 = C4
const double Const::c2 = (w0 + w1) / 2; // C2 = C3
const double Const::d1 = w1; // D1 = D3
const double Const::d2 = w0;