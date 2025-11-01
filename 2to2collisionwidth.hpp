#pragma once

#ifndef TWOTOTWOCOLLISIONWIDTH_HPP
#define TWOTOTWOCOLLISIONWIDTH_HPP

extern const int NDIM;
extern const int NCOMP;
extern const int NVEC;
extern const double EPSREL;
extern const double EPSABS;
extern const int VERBOSE;
extern const int LAST;
extern const int SEED;
extern const int MINEVAL;
extern const int MAXEVAL;
extern const int NSTART;
extern const int NINCREASE;
extern const int NBATCH;
extern const int GRIDNO;
extern const char *STATEFILE;
extern const int *SPIN;

#include "constants.hpp"
#include "numerics.hpp"
#include <algorithm>
#include <assert.h>
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

constexpr double _4_M_4        = 8 * M_PI * M_PI * M_PI; /* 8 pi^3 */
constexpr double _2_PI_FACTORS = pow(2 * M_PI, 4) / pow(2 * M_PI, 3 * 3);

constexpr double MAX_OMEGA_PLOL_EXP = log10(2);
using namespace boost::math::interpolators;

enum class ODDNESS
{
  EVEN,
  ODD,
};

template <typename T> int sgn(T val)
{
  return (T(0) < val) - (val < T(0));
}

namespace TwoToTwoCollisionWidth
{

/**
 * @brief vector addition
 */
template <typename T>
std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b)
{
  assert(a.size() == b.size());
  std::vector<T> result;
  result.reserve(a.size());

  std::transform(a.begin(),
                 a.end(),
                 b.begin(),
                 std::back_inserter(result),
                 std::plus<T>());
  return result;
}

/**
 * @brief vector subtraction
 */
template <typename T>
std::vector<T> operator-(const std::vector<T> &a, const std::vector<T> &b)
{
  assert(a.size() == b.size());

  std::vector<T> result;
  result.reserve(a.size());

  std::transform(a.begin(),
                 a.end(),
                 b.begin(),
                 std::back_inserter(result),
                 std::minus<T>());
  return result;
}

/**
 * @brief multiplication of vector with scalar
 */
template <typename T, typename T2>
std::vector<T> operator*(const T2 &a, const std::vector<T> &b)
{
  std::vector<T> result;
  result.reserve(b.size());

  std::transform(b.begin(),
                 b.end(),
                 std::back_inserter(result),
                 [&a](T i) { return a * i; });
  return result;
}

/**
 * @brief division of vector by scalar
 */
template <typename T, typename T2>
std::vector<T> operator/(const std::vector<T> &a, const T2 &b)
{
  std::vector<T> result;
  result.reserve(a.size());

  std::transform(a.begin(),
                 a.end(),
                 std::back_inserter(result),
                 [&b](T i) { return i / b; });
  return result;
}

/**
 * @brief dot product of two vectors
 */
template <typename T>
T operator*(const std::vector<T> &a, const std::vector<T> &b)
{
  assert(a.size() == b.size());

  std::vector<T> result;
  result.reserve(a.size());

  std::transform(a.begin(),
                 a.end(),
                 b.begin(),
                 std::back_inserter(result),
                 [](T i, T j) { return (i * j); });

  T result1 = std::accumulate(result.begin(), result.end(), 0.0);

  return result1;
}

/**
 * @brief multiplication of matrix with vector
 */
template <typename T>
std::vector<T> operator*(const std::vector<std::vector<T>> &a,
                         const std::vector<T> &b)
{
  if (a.size() != b.size())
    throw("Multiplication of matrix with vector cannot be done. Must have the "
          "same size.");

  std::vector<T> result;
  result.reserve(a.size());

  std::transform(a.begin(),
                 a.end(),
                 std::back_inserter(result),
                 [&](std::vector<T> i) { return (i * b); });

  return result;
}

template <typename F> class gsl_function_pp : public gsl_function
{
public:
  gsl_function_pp(const F &func) : _func(func)
  {
    function = &gsl_function_pp::invoke;
    params   = this;
  }

private:
  const F &_func;
  static double invoke(double x, void *params)
  {
    return static_cast<gsl_function_pp *>(params)->_func(x);
  }
};

struct Process
{

  std::vector<double> p1p2_, p1p3_;
  std::vector<double> p1_, p2_, p3_;
  double E1_, ET_, ED_, p12_, p13_;
  const int s1, s2, s3, s4;
  const double m1, m2, m3, m4, prefactor;
  double T;

  /**
   * @brief Calculate energy
   *
   * @param m mass
   * @param p momentum
   * @return double
   */
  static double Energy(const double &m, const std::vector<double> &p);

  /**
   * @brief Calculate the probability density
   *
   * @param m mass
   * @param p momentum
   * @param s -1 for boson, 1 for fermion
   * @return double
   */
  double
  Distribution(const double &m, const std::vector<double> &p, const int &s);

  /**
   * @brief \f$ L_1 \f$ function
   *
   * @param p momentum integrated
   * @param omega external energy
   * @param k external momentum
   * @return double
   */
  double L1(const double &p, const double &omega, const double &k);

  /**
   * @brief \f$ L_2 \f$ function
   *
   * @param p momentum integrated
   * @param omega external energy
   * @param k external momentum
   * @return double
   */
  double L2(const double &p, const double &omega, const double &k);

  /**
   * @brief Real part of \f$ T_1 \f$ function
   *
   * @param g Boson coupling constant
   * @param C \f$ C(R) \f$ of the fermion representation
   * @param omega incoming energy
   * @param k incoming momenta
   * @param mB boson mass
   * @param mF fermion mass
   * @return double
   */
  double ReT1(const double &g,
              const double &C,
              const double &omega,
              const double &k,
              const double &mB,
              const double &mF);

  /**
   * @brief Imaginary part of \f$ T_1 \f$ function
   *
   * @param g Boson coupling constant
   * @param C \f$ C(R) \f$ of the fermion representation
   * @param omega incoming energy
   * @param k incoming momenta
   * @param mB boson mass
   * @param mF fermion mass
   * @return double
   */
  double ImT1(const double &g,
              const double &C,
              const double &omega,
              const double &k,
              const double &mB,
              const double &mF);

  /**
   * @brief Real part of \f$ T_2 \f$ function
   *
   * @param g Boson coupling constant
   * @param C \f$ C(R) \f$ of the fermion representation
   * @param omega incoming energy
   * @param k incoming momenta
   * @param mB boson mass
   * @param mF fermion mass
   * @return double
   */
  double ReT2(const double &g,
              const double &C,
              const double &omega,
              const double &k,
              const double &mB,
              const double &mF);

  /**
   * @brief Imaginary part of \f$ T_2 \f$ function
   *
   * @param g Boson coupling constant
   * @param C \f$ C(R) \f$ of the fermion representation
   * @param omega incoming energy
   * @param k incoming momenta
   * @param mB boson mass
   * @param mF fermion mass
   * @return double
   */
  double ImT2(const double &g,
              const double &C,
              const double &omega,
              const double &k,
              const double &mB,
              const double &mF);

  /**
   * @brief Calculate r for a given rcoeff and shift. theta and phi are not used
   * in this routine
   *
   * @param r radius
   * @param delta_r derivative of inside delta function
   * @param shift shift to center the potato
   * @param theta angle with the x axis
   * @param phi angle with the z acis
   * @param rcoeff (p1 + p2) . director( p4 )
   */
  void calculate_r(double &r,
                   double &delta_r,
                   double shift,
                   double theta,
                   double phi,
                   const double &rcoeff);

  /**
   * @brief Integrates over the outgoing momentum for a given p1 and p2.
   * Integrates over the out going angles
   * @param E1 energy of particle 1
   * @param E2 energy of particle 2
   * @param p1 momentum of particle 1
   * @param p2 momentum of particle 2
   * @return double
   */
  double MonteCarloInt_s(const double &E1,
                         const double &E2,
                         const std::vector<double> p1,
                         const std::vector<double> p2);

  /**
   * @brief Integrand over theta and phi
   *
   * @param theta x axis angle
   * @param phi z axis angle
   * @return double
   */
  double integrand_theta_phi_s(const double &theta, const double &phi);

  /**
   * @brief Integrates over phi for a given theta
   *
   * @param theta
   * @return double
   */
  double integrate_phi_s(const double &theta);

  /**
   * @brief
   *
   * @param E1 Energy of particle 1
   * @param ET Energy of particle 1 + Energy of particle 2
   * @param p1 3-momentum of p1
   * @param p2 3-momentum of p2
   * @param p1p2 p1p2 = p1 + p2
   * @return double
   */
  double Integrate_s(const double &E1,
                     const double &ET,
                     const std::vector<double> &p1,
                     const std::vector<double> &p2,
                     const std::vector<double> &p1p2);

  /**
   * @brief Integrand over theta and phi
   *
   * @param r distance to z axis
   * @param phi z axis angle
   * @return double
   */
  double integrand_r_theta_t(const double &r, const double &theta);

  /**
   * @brief Integrates over phi for a given theta
   *
   * @param theta
   * @return double
   */
  double integrate_theta_t(const double &theta);
  /**
   * @brief Integrates over the outgoing momentum for a given p1 and p2.
   * Integrates over the out going angles
   * @param E1 energy of particle 1
   * @param E3 energy of particle 3
   * @param p1 momentum of particle 1
   * @param p3 momentum of particle 3
   * @return double
   */
  double MonteCarloInt_t(const double &E1,
                         const double &E3,
                         const std::vector<double> p1,
                         const std::vector<double> p3);

  /**
   * @brief
   *
   * @param E1 Energy of particle 1
   * @param ET Energy of particle 1 - Energy of particle 3
   * @param p1 3-momentum of p1
   * @param p3 3-momentum of p3
   * @param p1p3 p1p3 = p1 - p3
   * @return double
   */
  double Integrate_t(const double &E1,
                     const double &ED,
                     const std::vector<double> &p1,
                     const std::vector<double> &p3,
                     const std::vector<double> &p1p3);

  /**
   * @brief Constructor
   *
   * @param T_in Temperature
   * @param prefactor_in prefactor
   * @param s1_in = -1 boson, 1 for fermion. particle 1
   * @param s2_in = -1 boson, 1 for fermion. particle 2
   * @param s3_in = -1 boson, 1 for fermion. particle 3
   * @param s4_in = -1 boson, 1 for fermion. particle 4
   * @param m1_in mass of particle 1
   * @param m2_in mass of particle 2
   * @param m3_in mass of particle 3
   * @param m4_in mass of particle 4
   */
  explicit Process(const double &T_in,
                   const double &prefactor_in,
                   const int &s1_in,
                   const int &s2_in,
                   const int &s3_in,
                   const int &s4_in,
                   const double &m1_in,
                   const double &m2_in,
                   const double &m3_in,
                   const double &m4_in);

  /**
   * @brief propagator of channel_s
   *
   * @param p1 momentum of particle 1
   * @param p2 momentum of particle 2
   * @return double
   */
  virtual double PropagatorSquared_s(const std::vector<double> &p1,
                                     const std::vector<double> &p2) = 0;

  /**
   * @brief Amplitude of s channel
   *
   * @param p1 momentum of particle 1
   * @param p2 momentum of particle 2
   * @param p3 momentum of particle 3
   * @return double
   */
  virtual double AmplitudeSquared_s(const std::vector<double> &p1,
                                    const std::vector<double> &p2,
                                    const std::vector<double> &p3) = 0;

  /**
   * @brief propagator of channel_s
   *
   * @param p1 momentum of particle 1
   * @param p2 momentum of particle 2
   * @return double
   */
  virtual double PropagatorSquared_t(const std::vector<double> &p1,
                                     const std::vector<double> &p2) = 0;

  /**
   * @brief Amplitude of t channel
   *
   * @param p1 momentum of particle 1
   * @param p2 momentum of particle 2
   * @param p3 momentum of particle 3
   * @return double
   */
  virtual double AmplitudeSquared_t(const std::vector<double> &p1,
                                    const std::vector<double> &p2,
                                    const std::vector<double> &p3) = 0;

  static int Integrand_s(const int *ndim,
                         const cubareal xx[],
                         const int *ncomp,
                         cubareal ff[],
                         void *userdata);

  static int Integrand_t(const int *ndim,
                         const cubareal xx[],
                         const int *ncomp,
                         cubareal ff[],
                         void *userdata);
};
} // namespace TwoToTwoCollisionWidth
#endif
