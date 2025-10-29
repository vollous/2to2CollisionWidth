#pragma once

#include <boost/math/interpolators/bilinear_uniform.hpp>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdlib> // for std::getenv
#include <cuba.h>
#include <functional>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>

double Polylog(const int &s, const double &x);