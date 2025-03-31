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

#include <algorithm>
#include <assert.h>
#include <chrono>
#include <cmath>
#include <cuba.h>
#include <functional>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>

constexpr double _4_M_4        = 4 * M_PI * M_PI * M_PI * M_PI; /* 4 pi^4 */
constexpr double _2_PI_FACTORS = pow(2 * M_PI, 3) / pow(2 * M_PI, 3 * 4);

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

  std::vector<double> p1p2_;
  std::vector<double> p1_, p2_;
  double E1_, ET_, p12_;
  const int s1, s2, s3, s4;
  const double m1, m2, m3, m4, T;

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
  double MonteCarloInt(const double &E1,
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
  double integrand_theta_phi(const double &theta, const double &phi);

  static int Integrand(const int *ndim,
                       const cubareal xx[],
                       const int *ncomp,
                       cubareal ff[],
                       void *userdata);

  /**
   * @brief Integrates over phi for a given theta
   *
   * @param theta
   * @return double
   */
  double integrate_phi(const double &theta);

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
  double Integrate(const double &E1,
                   const double &ET,
                   const std::vector<double> &p1,
                   const std::vector<double> &p2,
                   const std::vector<double> &p1p2);

  /**
   * @brief Constructor
   *
   * @param T_in Temperature
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
                   const int &s1_in,
                   const int &s2_in,
                   const int &s3_in,
                   const int &s4_in,
                   const double &m1_in,
                   const double &m2_in,
                   const double &m3_in,
                   const double &m4_in);

  /**
   * @brief Amplitude
   *
   * @param p1 momentum of particle 1
   * @param p2 momentum of particle 2
   * @param p3 momentum of particle 3
   * @return double
   */
  virtual double AmplitudeSquared(const std::vector<double> &p1,
                                  const std::vector<double> &p2,
                                  const std::vector<double> &p3) = 0;
};
} // namespace TwoToTwoCollisionWidth
#endif
