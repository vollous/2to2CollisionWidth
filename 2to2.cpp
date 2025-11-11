#include "2to2collisionwidth.hpp"
#include "constants.hpp"
#include "matplotlibcpp.h"

const int NDIM        = 3;
const int NCOMP       = 1;
const int NVEC        = 1;
const double EPSREL   = 1e-2;
const double EPSABS   = 1e-12;
const int VERBOSE     = 3;
const int LAST        = 0;
const int SEED        = 0;
const int MINEVAL     = 0;
const int MAXEVAL     = 200000;
const int NSTART      = 300;
const int NINCREASE   = 300;
const int NBATCH      = 300;
const int GRIDNO      = 0;
const char *STATEFILE = NULL;
const int *SPIN       = NULL;

#define NNEW 1000
#define NMIN 2
#define FLATNESS 25.

#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0
using namespace boost::math::interpolators;
using namespace TwoToTwoCollisionWidth;
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;
namespace plt = matplotlibcpp;

struct tLgTOtRH_massless_helicity_thermal_masses : Process
{
  using Process::Process;                 // Import constructor
  const double mtinf = gs / sqrt(6.) * T; // Top thermal mass

  double AmplitudeSquared_s(const std::vector<double> &p1,
                            const std::vector<double> &p2,
                            const std::vector<double> &p3) override
  {
    return 0;
  }

  double AmplitudeSquared_t(const std::vector<double> &p1,
                            const std::vector<double> &p2,
                            const std::vector<double> &p3) override
  {

    const double t       = -2 * Energy(0, p1) * Energy(0, p3) + 2 * p1 * p3;
    const double p1dotp2 = p1 * p2;
    const double s       = 2 * Energy(0, p1) * Energy(0, p2) - 2 * p1dotp2;

    double res = (-2 * pow(el, 2) * pow(gs, 2) * pow(mt_pole, 2) * s * t) /
                 (pow(mW, 2) * pow(sW, 2)) / pow(-pow(mtinf, 2) + t, 2);

    if (res < 0)
    {
      std::cout << "Negative amplitude sqr? \n";
      return 0.;
    }
    return res;
  }
};

struct tLgTOtRH_massless_helicity_HTL : Process
{
  using Process::Process;           // Import constructor
  double mtinf = gs / sqrt(6.) * T; // Top thermal mass

  // Calculation of "a" and "b"

  inline double Rea(const double &m, const double &omega, const double &k)
  {
    const double logomegak = log(abs((omega + k) / (omega - k)));
    return pow(m / k, 2) * (1 - omega / (2. * k) * logomegak);
  }

  inline double Ima(const double &m, const double &omega, const double &k)
  {
    const double arglog = std::arg((omega + k) / (omega - k));
    return pow(m / k, 2) * (-omega / (2. * k) * arglog);
  }

  inline double Reb(const double &m, const double &omega, const double &k)
  {
    const double logomegak = log(abs((omega + k) / (omega - k)));
    return pow(m, 2) / k *
           (-omega / k + (pow(omega / k, 2) - 1) / 2. * logomegak);
  }

  inline double Imb(const double &m, const double &omega, const double &k)
  {
    const double arglog = std::arg((omega + k) / (omega - k));
    return pow(m, 2) / k * ((pow(omega / k, 2) - 1) / 2. * arglog);
  }

  double AmplitudeSquared_s(const std::vector<double> &p1,
                            const std::vector<double> &p2,
                            const std::vector<double> &p3) override
  {
    return 0;
  }

  double AmplitudeSquared_t(const std::vector<double> &p1,
                            const std::vector<double> &p2,
                            const std::vector<double> &p3) override
  {

    const double p1dotp2 = p1 * p2;
    const double s       = 2 * Energy(0, p1) * Energy(0, p2) - 2 * p1dotp2;
    const double t =
        -2 * Energy(0, p1) * Energy(0, p3) + 2 * p1 * p3; // omega^2-k^2

    // Self energy
    const std::complex<double> a(REa_t, IMa_t);
    const std::complex<double> b(REb_t, IMb_t);
    const double omega = Energy(0, p1) - Energy(0, p3);
    const double k     = Energy(0, p1 - p3);

    const double res =
        -(2 * pow(el, 2) * pow(gs, 2) * pow(mt_pole, 2) * s * t *
          std::norm((a + 1.) / (pow(b, 2) + 2. * (a + 1.) * b * omega +
                                pow(a + 1., 2) * t))) /
        (pow(mW, 2) * pow(sW, 2));

    if (res < 0)
    {
      std::cout << "Negative amplitude sqr? \n";
      return 0.;
    }
    return res;
  }
};

struct tLgTOtRH_massless_helicity_full : Process
{
  using Process::Process;                 // Import constructor
  const double mtinf = gs / sqrt(6.) * T; // Top thermal mass
  const double C     = 4. / 3.;

  std::function<double(double, double)> ReT1_rasterized;
  std::function<double(double, double)> ImT1_rasterized;
  std::function<double(double, double)> ReT2_rasterized;
  std::function<double(double, double)> ImT2_rasterized;

  inline std::function<double(double, double)>
  to_bilinear(const ODDNESS &oddness,
              const int &n1,
              const int &n2,
              const std::function<double(double)> &u,
              const std::function<double(double)> &v,
              const std::function<double(double, double)> &norm,
              const double &f_0_0,
              const std::function<double(double)> &f_v1_0,
              const double &f_inf_0,
              const std::function<double(double)> &f_0_v2,
              const std::function<double(double, double)> &f_v1_v2,
              const std::function<double(double)> &f_inf_v2,
              const double &f_0_inf,
              const std::function<double(double)> &f_v1_inf,
              const double &f_inf_inf)
  {
    if (u(v(M_PI)) / M_PI - 1 > 1e-6 or v(u(M_PI)) / M_PI - 1 > 1e-6)
      std::cout << "Error with u to v converter\n";

    if (n1 == n2)
      std::cout << "n1 = n2. Might be dangerous for some functions!";

    const double d1 = 1. / (n1 - 1.);
    const double d2 = 1. / (n2 - 1.);

    std::vector<double> vv;
    vv.push_back(f_0_0);
    for (double d_1 = d1; d_1 < 1 - d1 / 2.; d_1 += d1)
      vv.push_back(f_v1_0(v(d_1)));
    vv.push_back(f_inf_0);

    if (vv.size() != n1)
      std::cout << "Bad size at n1." << vv.size() << "\t" << n1 << "\n";

    for (double d_2 = d2; d_2 < 1 - d2 / 2.; d_2 += d2)
    {
      vv.push_back(f_0_v2(v(d_2)));
      for (double d_1 = d1; d_1 < 1 - d1 / 2.; d_1 += d1)
        vv.push_back(f_v1_v2(v(d_1), v(d_2)));
      vv.push_back(f_inf_v2(v(d_2)));
    }

    if (vv.size() != n1 * (n2 - 1))
      std::cout << "Bad size at n1.\t\t\t" << vv.size() << "\t" << n1 * (n2 - 1)
                << "\n";
    vv.push_back(f_0_inf);
    for (double d_1 = d1; d_1 < 1 - d1 / 2.; d_1 += d1)
    {
      // std::cout << "d_1\t" << d_1 << "\t" << v(d_1) << "\n";
      vv.push_back(f_v1_inf(v(d_1)));
    }

    vv.push_back(f_inf_inf);

    for (auto value : vv)
      if (isnan(value)) std::cout << "NaN value found\n";

    auto bu = bilinear_uniform(std::move(vv), n2, n1, d1, d2);
    std::function<double(double, double)> r;
    if (oddness == ODDNESS::ODD)
    {
      r = [=](const double &v1, const double &v2)
      { return sgn(v1) * bu(u(abs(v1)), u(v2)) / norm(abs(v1), v2); };
    }
    else
    {
      // even function first argument
      r = [=](const double &v1, const double &v2)
      { return bu(u(abs(v1)), u(v2)) / norm(abs(v1), v2); };
    }

    return r;
  }

  inline void Generate_Billinear_a_b(double mB = 0, double mF = 0)
  {
    // For testing purposes
    double oo = 10.1;
    double kk = 1.;

    // Hyper parameters
    const int n1       = 201;
    const int n2       = 202;
    const double s     = T;
    const double INFTY = 1.e5;
    const double ZERO  = 1.e-5;

    // var_tilde = s * var / (s * var + 1 )
    // var = var_tilde * (s - s var_tilde )
    std::function<double(double)> u = [=](const double &var)
    { return s * var / (1. + s * var); };
    std::function<double(double)> v = [=](const double &tilde)
    { return tilde / (1. - tilde) / s; };

    /*************************** ReT1 ***************************/

    std::function<double(double, double)> norm_ret1 =
        [=](const double &omega, const double &k)

    { return (omega - k + mF) / (omega + k + 1); };
    const double ret1_0_0                            = 0;
    const std::function<double(double)> ret1_omega_0 = [=](const double &omega)
    { return norm_ret1(omega, ZERO) * ReT1(gs, C, omega, ZERO, mB, mF); };

    const double ret1_inf_0 =
        norm_ret1(INFTY, ZERO) * ReT1(gs, C, INFTY, ZERO, mB, mF);

    const std::function<double(double)> ret1_0_k = [=](const double &k)
    { return norm_ret1(ZERO, k) * ReT1(gs, C, ZERO, k, mB, mF); };
    const std::function<double(double, double)> ret1_omega_k =
        [=](const double &omega, const double &k)
    { return norm_ret1(omega, k) * ReT1(gs, C, omega, k, mB, mF); };
    const std::function<double(double)> ret1_inf_k = [=](const double &k)
    { return norm_ret1(INFTY, k) * ReT1(gs, C, INFTY, k, mB, mF); };

    const double ret1_0_inf =
        norm_ret1(ZERO, INFTY) * ReT1(gs, C, ZERO, INFTY, mB, mF);
    const std::function<double(double)> ret1_omega_inf =
        [=](const double &omega)
    { return norm_ret1(omega, INFTY) * ReT1(gs, C, omega, INFTY, mB, mF); };
    const double ret1_inf_inf =
        norm_ret1(INFTY, INFTY + 1) * ReT1(gs, C, INFTY, INFTY + 1, mB, mF);

    ReT1_rasterized = to_bilinear(ODDNESS::ODD,
                                  n1,
                                  n2,
                                  u,
                                  v,
                                  norm_ret1,
                                  ret1_0_0,
                                  ret1_omega_0,
                                  ret1_inf_0,
                                  ret1_0_k,
                                  ret1_omega_k,
                                  ret1_inf_k,
                                  ret1_0_inf,
                                  ret1_omega_inf,
                                  ret1_inf_inf);

    std::cout << "ReT1_rasterized\t"
              << ReT1_rasterized(oo, kk) / ReT1(gs, C, oo, kk, mB, mF) - 1
              << "\n";

    /*************************** ImT1 ***************************/

    std::function<double(double, double)> norm_imt1 =
        [=](const double &omega, const double &k)
    { return k * (omega - k + mF) / (omega + k + 1); };

    const double imt1_0_0                            = 0;
    const std::function<double(double)> imt1_omega_0 = [=](const double &omega)
    { return norm_imt1(omega, ZERO) * ImT1(gs, C, omega, ZERO, mB, mF); };
    const double imt1_inf_0 = 0;

    const std::function<double(double)> imt1_0_k = [=](const double &k)
    {
      return (
          (pow(k, 2) * T *
           (-(sqrt(pow(k, 2) + 4 * pow(mB, 2)) *
              log(1 - exp(-0.5 * sqrt(pow(k, 2) + 4 * pow(mB, 2)) / T))) +
            sqrt(pow(k, 2) + 4 * pow(mF, 2)) *
                log(1 + exp(-0.5 * sqrt(pow(k, 2) + 4 * pow(mF, 2)) / T)) +
            2 * T *
                Polylog(2, exp(-0.5 * sqrt(pow(k, 2) + 4 * pow(mB, 2)) / T)) -
            2 * T *
                Polylog(2,
                        -exp(-0.5 * sqrt(pow(k, 2) + 4 * pow(mF, 2)) / T)))) /
          (1 + k));
    };
    const std::function<double(double, double)> imt1_omega_k =
        [=](const double &omega, const double &k)
    { return norm_imt1(omega, k) * ImT1(gs, C, omega, k, mB, mF); };
    const std::function<double(double)> imt1_inf_k = [=](const double &k)
    { return 0; };

    const double imt1_0_inf = 0;
    const std::function<double(double)> imt1_omega_inf =
        [=](const double &omega) { return 0; };
    const double imt1_inf_inf = 0;

    ImT1_rasterized = to_bilinear(ODDNESS::EVEN,
                                  n1,
                                  n2,
                                  u,
                                  v,
                                  norm_imt1,
                                  imt1_0_0,
                                  imt1_omega_0,
                                  imt1_inf_0,
                                  imt1_0_k,
                                  imt1_omega_k,
                                  imt1_inf_k,
                                  imt1_0_inf,
                                  imt1_omega_inf,
                                  imt1_inf_inf);

    std::cout << "ImT1_rasterized\t"
              << ImT1_rasterized(oo, kk) / ImT1(gs, C, oo, kk, mB, mF) - 1
              << "\n";

    /*************************** ReT2 ***************************/

    std::function<double(double, double)> norm_ret2 =
        [=](const double &omega, const double &k) { return 1; };
    const double ret2_0_0                            = 0;
    const std::function<double(double)> ret2_omega_0 = [=](const double &omega)
    { return norm_ret2(omega, ZERO) * ReT2(gs, C, omega, ZERO, mB, mF); };

    const double ret2_inf_0 =
        norm_ret2(INFTY, ZERO) * ReT2(gs, C, INFTY, ZERO, mB, mF);

    const std::function<double(double)> ret2_0_k = [=](const double &k)
    { return norm_ret2(ZERO, k) * ReT2(gs, C, ZERO, k, mB, mF); };
    const std::function<double(double, double)> ret2_omega_k =
        [=](const double &omega, const double &k)
    { return norm_ret2(omega, k) * ReT2(gs, C, omega, k, mB, mF); };
    const std::function<double(double)> ret2_inf_k = [=](const double &k)
    { return norm_ret2(INFTY, k) * ReT2(gs, C, INFTY, k, mB, mF); };

    const double ret2_0_inf =
        norm_ret2(ZERO, INFTY) * ReT2(gs, C, ZERO, INFTY, mB, mF);
    const std::function<double(double)> ret2_omega_inf =
        [=](const double &omega)
    { return norm_ret2(omega, INFTY) * ReT2(gs, C, omega, INFTY, mB, mF); };
    const double ret2_inf_inf = 0;

    ReT2_rasterized = to_bilinear(ODDNESS::EVEN,
                                  n1,
                                  n2,
                                  u,
                                  v,
                                  norm_ret2,
                                  ret2_0_0,
                                  ret2_omega_0,
                                  ret2_inf_0,
                                  ret2_0_k,
                                  ret2_omega_k,
                                  ret2_inf_k,
                                  ret2_0_inf,
                                  ret2_omega_inf,
                                  ret2_inf_inf);

    std::cout << "ReT2_rasterized\t"
              << ReT2_rasterized(oo, kk) / ReT2(gs, C, oo, kk, mB, mF) - 1
              << "\n";

    /*************************** ImT2 ***************************/

    std::function<double(double, double)> norm_imt2 =
        [=](const double &omega, const double &k)

    { return (omega - k + mF) / (omega + k + 1); };
    const double imt2_0_0                            = 0;
    const std::function<double(double)> imt2_omega_0 = [=](const double &omega)
    { return norm_imt2(omega, ZERO) * ImT2(gs, C, omega, ZERO, mB, mF); };

    const double imt2_inf_0 =
        norm_imt2(INFTY, ZERO) * ImT2(gs, C, INFTY, ZERO, mB, mF);

    const std::function<double(double)> imt2_0_k = [=](const double &k)
    { return norm_imt2(ZERO, k) * ImT2(gs, C, ZERO, k, mB, mF); };
    const std::function<double(double, double)> imt2_omega_k =
        [=](const double &omega, const double &k)
    { return norm_imt2(omega, k) * ImT2(gs, C, omega, k, mB, mF); };
    const std::function<double(double)> imt2_inf_k = [=](const double &k)
    { return norm_imt2(INFTY, k) * ImT2(gs, C, INFTY, k, mB, mF); };

    const double imt2_0_inf =
        norm_imt2(ZERO, INFTY) * ImT2(gs, C, ZERO, INFTY, mB, mF);
    const std::function<double(double)> imt2_omega_inf =
        [=](const double &omega)
    { return norm_imt2(omega, INFTY) * ImT2(gs, C, omega, INFTY, mB, mF); };
    const double imt2_inf_inf =
        norm_imt2(INFTY, INFTY + 1) * ImT2(gs, C, INFTY, INFTY + 1, mB, mF);

    ImT2_rasterized = to_bilinear(ODDNESS::ODD,
                                  n1,
                                  n2,
                                  u,
                                  v,
                                  norm_imt2,
                                  imt2_0_0,
                                  imt2_omega_0,
                                  imt2_inf_0,
                                  imt2_0_k,
                                  imt2_omega_k,
                                  imt2_inf_k,
                                  imt2_0_inf,
                                  imt2_omega_inf,
                                  imt2_inf_inf);

    std::cout << "ImT2_rasterized\t"
              << ImT2_rasterized(oo, kk) / ImT2(gs, C, oo, kk, mB, mF) - 1
              << "\n";

    // rasterized_rasterized exit(0);
  }

  tLgTOtRH_massless_helicity_full(const double &T_in,
                                  const double &prefactor_in,
                                  const int &s1_in,
                                  const int &s2_in,
                                  const int &s3_in,
                                  const int &s4_in,
                                  const double &m1_in,
                                  const double &m2_in,
                                  const double &m3_in,
                                  const double &m4_in)
      : Process(T_in,
                prefactor_in,
                s1_in,
                s2_in,
                s3_in,
                s4_in,
                m1_in,
                m2_in,
                m3_in,
                m4_in) {
        // Generate_Billinear_a_b();
      };

  // Calculation of "a" and "b"

  void Calculate_a_b(const double &omega,
                     const double &k,
                     double &REa,
                     double &IMa,
                     double &REb,
                     double &IMb) override
  {
    const double rT1 = ReT1(gs, C, omega, k, 0, 0);
    const double rT2 = ReT2(gs, C, omega, k, 0, 0);
    const double iT1 = ImT1(gs, C, omega, k, 0, 0);
    const double iT2 = ImT2(gs, C, omega, k, 0, 0);

    /*const double rT1 = ReT1_rasterized(omega, k);
    const double rT2 = ReT2_rasterized(omega, k);
    const double iT1 = ImT1_rasterized(omega, k);
    const double iT2 = ImT2_rasterized(omega, k);*/

    REa = rT2 / pow(k, 2) - (omega * rT1) / pow(k, 2);
    IMa = iT2 / pow(k, 2) - (omega * iT1) / pow(k, 2);
    REb = -((omega * rT2) / pow(k, 2)) -
          ((pow(k, 2) - pow(omega, 2)) * rT1) / pow(k, 2);
    IMb = -((omega * iT2) / pow(k, 2)) -
          ((pow(k, 2) - pow(omega, 2)) * iT1) / pow(k, 2);
    return;
  }

  double AmplitudeSquared_s(const std::vector<double> &p1,
                            const std::vector<double> &p2,
                            const std::vector<double> &p3) override
  {
    return 0;
  }

  double AmplitudeSquared_t(const std::vector<double> &p1,
                            const std::vector<double> &p2,
                            const std::vector<double> &p3) override
  {

    const double p1dotp2 = p1 * p2;
    const double s       = 2 * Energy(0, p1) * Energy(0, p2) - 2 * p1dotp2;
    const double t       = -2 * Energy(0, p1) * Energy(0, p3) + 2 * p1 * p3;

    const std::complex<double> a(REa_t, IMa_t);
    const std::complex<double> b(REb_t, IMb_t);
    const double omega = Energy(0, p1) - Energy(0, p3);
    const double k     = Energy(0, p1 - p3);

    const double res =
        -(2 * pow(el, 2) * pow(gs, 2) * pow(mt_pole, 2) * s * t *
          std::norm((a + 1.) / (pow(b, 2) + 2. * (a + 1.) * b * omega +
                                pow(a + 1., 2) * t))) /
        (pow(mW, 2) * pow(sW, 2));

    if (res < 0)
    {
      std::cout << "Negative amplitude sqr? \n";
      return 0.;
    }
    return res;
  }
};

/********************* tL g -> tR g *********************/

struct tLgTOtRg_massless_helicity : Process
{
  using Process::Process;           // Import constructor
  double mtinf = gs / sqrt(6.) * T; // Top thermal mass

  double AmplitudeSquared_s(const std::vector<double> &p1,
                            const std::vector<double> &p2,
                            const std::vector<double> &p3) override
  {
    return 0.;
  }

  double AmplitudeSquared_t(const std::vector<double> &p1,
                            const std::vector<double> &p2,
                            const std::vector<double> &p3) override
  {

    const double t       = -2 * Energy(0, p1) * Energy(0, p3) + 2 * p1 * p3;
    const double p1dotp2 = p1 * p2;
    const double s       = 2 * Energy(0, p1) * Energy(0, p2) - 2 * p1dotp2;

    double res = (64 * pow(gs, 4) * (s + t)) / (3. * pow(pow(mtinf, 2) - s, 2));

    if (res < 0)
    {
      std::cout << "Negative amplitude sqr? \n";
      return 0.;
    }
    return res;
  }
};

struct identity : Process
{
  using Process::Process; // Import constructor

  double AmplitudeSquared_s(const std::vector<double> &p1,
                            const std::vector<double> &p2,
                            const std::vector<double> &p3) override
  {
    return 0.;
  }

  double AmplitudeSquared_t(const std::vector<double> &p1,
                            const std::vector<double> &p2,
                            const std::vector<double> &p3) override
  {
    return 1.;
  }
};

void warm_up_vegas(integrand_t integrand,
                   Process &proc,
                   int points,
                   int iterations,
                   int gridno,
                   cubareal integral[],
                   cubareal error[],
                   cubareal prob[])
{
  int last_only = 1, smoothing = 1;
  int comp, nregions, neval, fail;
  int maxpoints  = iterations * points;
  int cuba_flags = VERBOSE + last_only * 4 + smoothing * 8;

  Vegas(NDIM,
        NCOMP,
        integrand,
        &proc,
        NVEC,
        EPSREL,
        EPSABS,
        cuba_flags,
        SEED,
        MINEVAL,
        maxpoints,
        points,
        NINCREASE,
        NBATCH,
        gridno,
        STATEFILE,
        &SPIN,
        &neval,
        &fail,
        integral,
        error,
        prob);
}

void gridded_vegas(integrand_t integrand,
                   Process &proc,
                   int points,
                   int iterations,
                   int gridno,
                   cubareal integral[],
                   cubareal error[],
                   cubareal prob[])
{

  int last_only = 0, smoothing = 1;
  int comp, nregions, neval, fail;
  int maxpoints  = iterations * points;
  int cuba_flags = VERBOSE + last_only * 4 + smoothing * 8;
  // char *state = "my_state";
  Vegas(NDIM,
        NCOMP,
        integrand,
        &proc,
        NVEC,
        EPSREL,
        EPSABS,
        cuba_flags,
        SEED,
        MINEVAL,
        maxpoints,
        points,
        NINCREASE,
        NBATCH,
        gridno,
        STATEFILE,
        &SPIN,
        &neval,
        &fail,
        integral,
        error,
        prob);
}
void test_ReT1()
{
  plt::cla();
  double T = 1.;
  double v = 10;
  double m1, m2, m3, m4;

  m1 = 0; // mt
  m2 = 0; // mg
  m3 = 0; // ms
  m4 = 0; // mt

  int s1 = 1;
  int s2 = -1;
  int s3 = 1;
  int s4 = -1;

  tLgTOtRH_massless_helicity_HTL proc(T, T * T, s1, s2, s3, s4, m1, m2, m3, m4);

  for (double k : {0.5, 1.0, 1.5})
  {
    std::vector<double> x, y;
    auto tt1       = high_resolution_clock::now();
    double maxtime = -1;
    for (double i = -6; i <= MAX_OMEGA_PLOL_EXP; i += 0.0001)
    {
      double omega = pow(10, i);
      x.push_back(omega);
      auto t1                                = high_resolution_clock::now();
      double val                             = proc.ReT1(1, 1, omega, k, 0, 0);
      auto t2                                = high_resolution_clock::now();
      duration<double, std::milli> ms_double = t2 - t1;
      maxtime = std::max(maxtime, ms_double.count());
      // std::cout << "Re(T1, k = " << k << ", omega = " << omega
      //           << ") took around " << ms_double.count() / 1000. << "s\n";
      y.push_back(val);
    }
    auto tt2                               = high_resolution_clock::now();
    duration<double, std::milli> ms_double = tt2 - tt1;
    std::cout << "Re(T1, k = " << k << ") took around " << ms_double.count()
              << " ms.\t Max evaluation = " << maxtime << " ms.\n";
    plt::plot(x, y, {{"label", "k=" + std::to_string(k)}});
  }
  plt::xlabel("$\\omega$");
  plt::ylabel("$\\Re(T_1)$");
  plt::legend();
  plt::tight_layout();
  plt::save("../plots/ret1.pdf");
}

void test_ImT1()
{
  plt::cla();
  double T = 1.;
  double v = 10;
  double m1, m2, m3, m4;

  m1 = 0; // mt
  m2 = 0; // mg
  m3 = 0; // ms
  m4 = 0; // mt

  int s1 = 1;
  int s2 = -1;
  int s3 = 1;
  int s4 = -1;

  tLgTOtRH_massless_helicity_HTL proc(T, T * T, s1, s2, s3, s4, m1, m2, m3, m4);

  for (double k : {0.5, 1.0, 1.5})
  {
    std::vector<double> x, y;
    auto tt1       = high_resolution_clock::now();
    double maxtime = -1;
    for (double i = -6; i <= MAX_OMEGA_PLOL_EXP; i += 0.00011)
    {
      double omega = pow(10, i);
      x.push_back(omega);
      auto t1                                = high_resolution_clock::now();
      double val                             = proc.ImT1(1, 1, omega, k, 0, 0);
      auto t2                                = high_resolution_clock::now();
      duration<double, std::milli> ms_double = t2 - t1;
      maxtime = std::max(maxtime, ms_double.count());
      // std::cout << "Im(T1, k = " << k << ", omega = " << omega
      //           << ") took around " << ms_double.count() / 1000. << "s\n";
      y.push_back(val);
    }
    auto tt2                               = high_resolution_clock::now();
    duration<double, std::milli> ms_double = tt2 - tt1;
    std::cout << "Im(T1, k = " << k << ") took around " << ms_double.count()
              << " ms.\t Max evaluation = " << maxtime << " ms.\n";
    plt::plot(x, y, {{"label", "k=" + std::to_string(k)}});
  }
  plt::xlabel("$\\omega$");
  plt::ylabel("$\\Im(T_1)$");
  plt::legend();
  plt::tight_layout();
  plt::save("../plots/imt1.pdf");
}

void test_ReT2()
{
  plt::cla();
  double T = 1.;
  double v = 10;
  double m1, m2, m3, m4;

  m1 = 0; // mt
  m2 = 0; // mg
  m3 = 0; // ms
  m4 = 0; // mt

  int s1 = 1;
  int s2 = -1;
  int s3 = 1;
  int s4 = -1;

  tLgTOtRH_massless_helicity_HTL proc(T, T * T, s1, s2, s3, s4, m1, m2, m3, m4);

  for (double k : {0.5, 1.0, 1.5})
  {
    std::vector<double> x, y;
    auto tt1       = high_resolution_clock::now();
    double maxtime = -1;
    for (double i = -6; i <= MAX_OMEGA_PLOL_EXP; i += 0.0001)
    {
      double omega = pow(10, i);
      x.push_back(omega);
      auto t1                                = high_resolution_clock::now();
      double val                             = proc.ReT2(1, 1, omega, k, 0, 0);
      auto t2                                = high_resolution_clock::now();
      duration<double, std::milli> ms_double = t2 - t1;
      maxtime = std::max(maxtime, ms_double.count());
      // std::cout << "Re(T2, k = " << k << ", omega = " << omega
      //           << ") took around " << ms_double.count() / 1000. << "s\n";
      y.push_back(val);
    }
    auto tt2                               = high_resolution_clock::now();
    duration<double, std::milli> ms_double = tt2 - tt1;
    std::cout << "Re(T2, k = " << k << ") took around " << ms_double.count()
              << " ms.\t Max evaluation = " << maxtime << " ms.\n";
    plt::plot(x, y, {{"label", "k=" + std::to_string(k)}});
  }
  plt::xlabel("$\\omega$");
  plt::ylabel("$\\Re(T_2)$");
  plt::legend();
  plt::tight_layout();
  plt::save("../plots/ret2.pdf");
}

void test_ImT2()
{
  plt::cla();
  double T = 1.;
  double v = 10;
  double m1, m2, m3, m4;

  m1 = 0; // mt
  m2 = 0; // mg
  m3 = 0; // ms
  m4 = 0; // mt

  int s1 = 1;
  int s2 = -1;
  int s3 = 1;
  int s4 = -1;

  tLgTOtRH_massless_helicity_HTL proc(T, T * T, s1, s2, s3, s4, m1, m2, m3, m4);

  for (double k : {0.5, 1.0, 1.5})
  {
    std::vector<double> x, y;
    auto tt1       = high_resolution_clock::now();
    double maxtime = -1;
    for (double i = -6; i <= MAX_OMEGA_PLOL_EXP; i += 0.0001)
    {
      double omega = pow(10, i);
      x.push_back(omega);
      auto t1                                = high_resolution_clock::now();
      double val                             = proc.ImT2(1, 1, omega, k, 0, 0);
      auto t2                                = high_resolution_clock::now();
      duration<double, std::milli> ms_double = t2 - t1;
      maxtime = std::max(maxtime, ms_double.count());
      // std::cout << "Im(T2, k = " << k << ", omega = " << omega
      //           << ") took around " << ms_double.count() / 1000. << "s\n";
      y.push_back(val);
    }
    auto tt2                               = high_resolution_clock::now();
    duration<double, std::milli> ms_double = tt2 - tt1;
    std::cout << "Im(T2, k = " << k << ") took around " << ms_double.count()
              << " ms.\t Max evaluation = " << maxtime << " ms.\n";
    plt::plot(x, y, {{"label", "k=" + std::to_string(k)}});
  }
  plt::xlabel("$\\omega$");
  plt::ylabel("$\\Im(T_2)$");
  plt::legend();
  plt::tight_layout();
  plt::save("../plots/imt2.pdf");
}

void testing()
{

  std::cout
      << "\n\t\t\t ------------------ TESTING MODE ------------------\n\n";

  test_ReT1();
  std::cout << "\n";
  test_ImT1(); // No problem
  std::cout << "\n";
  test_ReT2();
  std::cout << "\n";
  test_ImT2(); // No problem
  std::cout << "\n";

  exit(0); // Exit early
}

int main()
{
  if (std::getenv("TESTING")) testing();

  /*int ncores = 5, pcores = 1e3;
  cubacores(&ncores, &pcores);*/

  double T = 10.;
  double v = 10;
  double m1, m2, m3, m4;

  std::cout << "Temperature is T = " << T << "\n";

  double mt = sqrt(pow(yt * v / sqrt(2.), 2) + pow(gs / sqrt(6.) * T, 2));
  double mg = gs * sqrt(3. / 3. + 3. / 6.) * T;
  double ms = 100;
  std::cout << "\n\nTop mass = " << mt << " GeV\n";
  std::cout << "Scalar mass = " << ms << " GeV\n";
  std::cout << "Gluon mass = " << mg << " GeV\n\n\n";

  int s1 = 1;
  int s2 = -1;
  int s3 = 1;
  int s4 = -1;

  /****************************************************************************/
  /******************************* tL g -> tR h *******************************/
  /****************************************************************************/

  /*
    t_L                     h
       \                  /
        \                /
         \              /
           ------------
         &              \
        &                \
       &                  \
      g                    t_R
*/

  // t-channel with massless approximationsi. // care very much about helicities
  // Gamma_y = 0.00571382 +- 3.30013e-05
  // Time : 40.739 s
  // Gamma_y = 0.056501 +- 0.000246043
  // Time redux : 30.132 s
  m1 = 0; // mt
  m2 = 0; // mg
  m3 = 0; // ms
  m4 = 0; // mt
  // tLgTOtRH_massless_helicity_thermal_masses proc(
  //     T, T * T, s1, s2, s3, s4, m1, m2, m3, m4);

  //  t-channel with massless HTL. // care very much about helicities
  //  Gamma_y = = 0.00422475 +- 1.73651e-05
  //  Time :17.969  s
  m1 = 0; // mt
  m2 = 0; // mg
  m3 = 0; // ms
  m4 = 0; // mt
  // tLgTOtRH_massless_helicity_HTL proc(T, T * T, s1, s2, s3, s4, m1, m2, m3,
  // m4);

  // std::vector<double> p1 = {-10, 5, 1};
  // std::vector<double> p3 = {10, 5, -2};
  //  std::cout << "->"
  //            << proc.MonteCarloInt_t(
  //                   proc.Energy(m1, p1), proc.Energy(m3, p3), p1, p3)
  //            << "\n";
  //   exit(0);

  //  t-channel with massless full propagator. // care very much about
  //  helicities Gamma_y = Time : Gamma_y = Time redux :
  m1 = 0; // mt
  m2 = 0; // mg
  m3 = 0; // ms
  m4 = 0; // mt
  tLgTOtRH_massless_helicity_full proc(
      T, T * T, s1, s2, s3, s4, m1, m2, m3, m4);

  /****************************** useless *******************************/

  // identity proc(100, 1, 0, 0, 0, 0, 0, 0, 0, 0);
  // identity proc(100, 1, s1, s2, s3, s4, 0, 0, 0, 0);

  int comp, nregions, neval, fail;
  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

  auto start = std::chrono::system_clock::now();
  int mode   = 1;
  if (mode == 0)
  {
    warm_up_vegas(proc.Integrand_t, proc, 1e3, 20, -1, integral, error, prob);
    gridded_vegas(proc.Integrand_t, proc, 1e4, 10, 1, integral, error, prob);
  }
  else
  {
    Vegas(NDIM,
          NCOMP,
          proc.Integrand_t,
          &proc,
          NVEC,
          EPSREL,
          EPSABS,
          VERBOSE,
          SEED,
          MINEVAL,
          MAXEVAL,
          NSTART,
          NINCREASE,
          NBATCH,
          GRIDNO,
          STATEFILE,
          &SPIN,
          &neval,
          &fail,
          integral,
          error,
          prob);
  }

  std::cout << "Integral\t" << integral[0] << "\n";
  std::cout << "Error \t" << error[0] << "\n";
  std::cout << "Relative error \t" << error[0] / integral[0] << "\n";
  std::cout << "prob\t" << prob[0] << "\n";

  const double vw    = 0.95;
  const double gamma = 1; /*/ sqrt(1 - pow(vw, 2));*/
  const double N1    = gamma * 2 * pow(M_PI, 3) * pow(T, 2) / 3;
  std::cout << "\nT = " << T
            << " | Gamma_y = " << integral[0] / (N1 * proc.prefactor) << " +- "
            << error[0] / (N1 * proc.prefactor) << "\n\n";

  // Your Code to Execute //
  auto end = std::chrono::system_clock::now();
  std::cout << "It took\t"
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     start)
                       .count() /
                   1000.
            << "\ts" << std::endl;
  return 0;
}
