#pragma once
#include <cmath>

constexpr const double alphaS  = 0.12;
constexpr const double ew      = 246.22;
constexpr const double el      = 0.332;
constexpr const double sW      = sqrt(0.223);
constexpr const double gs      = sqrt(4 * M_PI * alphaS);
constexpr const double mW      = 80.;
constexpr const double mt_pole = 172.5;

constexpr const double yt = mt_pole / (ew / sqrt(2.));