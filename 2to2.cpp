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
const int NSTART      = 2000;
const int NINCREASE   = 2000;
const int NBATCH      = 2000;
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
using namespace TwoToTwoCollisionWidth;
namespace plt = matplotlibcpp;
struct tLgTOtRH : Process
{
  using Process::Process; // Import constructor

  double AmplitudeSquared(const std::vector<double> &p1,
                          const std::vector<double> &p2,
                          const std::vector<double> &p3) override
  {

    const double t =
        pow(m1, 2) + pow(m3, 2) - 2 * E1_ * Energy(m3, p3) + 2 * p1_ * p3;
    const double p1dotp2 = p1 * p2;
    const double s = m1 * m1 + m2 * m2 + 2 * E1_ * Energy(m2, p2) - 2 * p1dotp2;
    const double mt = m1;
    const double mg = m2;
    const double ms = m3;
    return (
        (2 * pow(gs, 2) *
         (2 * pow(mg, 2) *
              (pow(mg, 4) * pow(s, 2) +
               pow(mg, 2) * (pow(ms, 2) * pow(s - t, 2) + pow(s, 2) * t) -
               s * t *
                   (2 * pow(ms, 4) - 2 * pow(ms, 2) * (s + t) +
                    pow(s + t, 2))) +
          pow(mt, 2) *
              (-4 * pow(mg, 6) * s - s * t * (pow(s, 2) + pow(t, 2)) -
               2 * pow(mg, 4) * (4 * pow(s, 2) - 2 * s * t + pow(t, 2)) +
               pow(mg, 2) * (4 * pow(ms, 4) * (s + t) + 8 * s * t * (s + t) -
                             5 * pow(ms, 2) * (pow(s, 2) + pow(t, 2))) +
               pow(mt, 2) *
                   (2 * pow(mg, 6) + pow(mg, 4) * (8 * s - 2 * t) +
                    pow(s + t, 3) +
                    pow(mt, 2) *
                        (-2 * pow(mg, 4) + 6 * pow(mg, 2) * pow(ms, 2) -
                         2 * pow(mt, 4) + 4 * pow(mt, 2) * (s + t) -
                         3 * pow(s + t, 2)) -
                    2 * pow(mg, 2) *
                        (2 * pow(ms, 4) + pow(ms, 2) * (s + t) -
                         2 * (pow(s, 2) - 4 * s * t + pow(t, 2)))))) *
         pow(yt, 2)) /
        (pow(mg, 2) * pow(-pow(mt, 2) + s, 2) * pow(-pow(mt, 2) + t, 2)));
  }
};

struct tLgTOtRH_massless_gluon : Process
{
  using Process::Process; // Import constructor

  double AmplitudeSquared(const std::vector<double> &p1,
                          const std::vector<double> &p2,
                          const std::vector<double> &p3) override
  {

    const double t =
        pow(m1, 2) + pow(m3, 2) - 2 * E1_ * Energy(m3, p3) + 2 * p1_ * p3;
    const double p1dotp2 = p1 * p2;
    const double s = m1 * m1 + m2 * m2 + 2 * E1_ * Energy(m2, p2) - 2 * p1dotp2;
    const double mt = m1;
    const double ms = m3;
    return (
        (4 * pow(gs, 2) *
         (2 * pow(ms, 2) * pow(mt, 6) +
          pow(mt, 4) * (-2 * pow(ms, 4) + (s - 3 * t) * (3 * s - t)) -
          s * t * (2 * pow(ms, 4) - 2 * pow(ms, 2) * (s + t) + pow(s + t, 2)) +
          pow(mt, 2) * (2 * pow(ms, 4) * (s + t) + 4 * s * t * (s + t) -
                        3 * pow(ms, 2) * (pow(s, 2) + pow(t, 2)))) *
         pow(yt, 2)) /
        pow(pow(mt, 4) + s * t - pow(mt, 2) * (s + t), 2));
  }
};

struct tLgTOtRH_massless : Process
{
  using Process::Process;           // Import constructor
  double mtinf = gs / sqrt(6.) * T; // Top thermal mass
  double AmplitudeSquared(const std::vector<double> &p1,
                          const std::vector<double> &p2,
                          const std::vector<double> &p3) override
  {

    const double t       = -2 * Energy(0, p1) * Energy(0, p3) + 2 * p1 * p3;
    const double p1dotp2 = p1 * p2;
    const double s       = 2 * Energy(0, p1) * Energy(0, p2) - 2 * p1dotp2;

    double res = (-4 * s * t * el * el * gs * gs * mt_pole * mt_pole) /
                 (3. * pow(t - mtinf * mtinf, 2) * mW * mW * sW * sW);

    if (res < 0)
    {
      std::cout << "Negative amplitude sqr? \n";
      return 0.;
    }
    return res;
  }
};

struct tLgTOtRH_massless_helicity_thermal_masses : Process
{
  using Process::Process;                 // Import constructor
  const double mtinf = gs / sqrt(6.) * T; // Top thermal mass
  double AmplitudeSquared(const std::vector<double> &p1,
                          const std::vector<double> &p2,
                          const std::vector<double> &p3) override
  {

    const double t       = -2 * Energy(0, p1) * Energy(0, p3) + 2 * p1 * p3;
    const double p1dotp2 = p1 * p2;
    const double s       = 2 * Energy(0, p1) * Energy(0, p2) - 2 * p1dotp2;

    double res = (-2 * pow(el, 2) * pow(gs, 2) * pow(mt_pole, 2) * s * t) /
                 (pow(mW, 2) * pow(sW, 2) * pow(-pow(mtinf, 2) + t, 2));

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
    return pow(m, 2) / k * (1 - omega / (2. * k) * logomegak);
  }

  inline double Ima(const double &m, const double &omega, const double &k)
  {
    const double arglog = std::arg((omega + k) / (omega - k));
    return pow(m, 2) / k * (-omega / (2. * k) * arglog);
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

  double AmplitudeSquared(const std::vector<double> &p1,
                          const std::vector<double> &p2,
                          const std::vector<double> &p3) override
  {

    const double p1dotp2 = p1 * p2;
    const double s       = 2 * Energy(0, p1) * Energy(0, p2) - 2 * p1dotp2;

    // Self energy
    const double omega = Energy(0, p1) - Energy(0, p3);
    const double k     = Energy(0, p1 - p3);
    const double t =
        -2 * Energy(0, p1) * Energy(0, p3) + 2 * p1 * p3; // omega^2-k^2

    const double REa = Rea(mtinf, omega, k);
    const double IMa = Ima(mtinf, omega, k);
    const double REb = Reb(mtinf, omega, k);
    const double IMb = Imb(mtinf, omega, k);

    double res =
        ((2 * pow(el, 2) * pow(gs, 2) * pow(mt_pole, 2) *
          (pow(IMa, 2) + pow(1 + REa, 2))) /
         (pow(mW, 2) * pow(sW, 2) *
          (pow(IMb, 4) + 4 * IMa * pow(IMb, 3) * omega +
           4 * omega * (1 + REa) * pow(REb, 3) + pow(REb, 4) +
           4 * omega * (1 + REa) * (pow(IMa, 2) + pow(1 + REa, 2)) * REb * t +
           pow(pow(IMa, 2) + pow(1 + REa, 2), 2) * pow(t, 2) +
           2 * pow(IMb, 2) *
               (2 * pow(omega, 2) * (pow(IMa, 2) + pow(1 + REa, 2)) +
                2 * omega * (1 + REa) * REb + pow(REb, 2) +
                (pow(IMa, 2) - pow(1 + REa, 2)) * t) +
           2 * pow(REb, 2) *
               (2 * pow(omega, 2) * (pow(IMa, 2) + pow(1 + REa, 2)) +
                (-pow(IMa, 2) + pow(1 + REa, 2)) * t) +
           4 * IMa * IMb *
               (omega * pow(REb, 2) +
                omega * (pow(IMa, 2) + pow(1 + REa, 2)) * t +
                2 * (1 + REa) * REb * t))));

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
  double AmplitudeSquared(const std::vector<double> &p1,
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

  double AmplitudeSquared(const std::vector<double> &p1,
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

  tLgTOtRH_massless_gluon proc(T, T * T, s1, s2, s3, s4, m1, m2, m3, m4);

  for (double k : {0.5, 1.0, 1.5})
  {
    std::vector<double> x, y;

    for (double omega = 0; omega < 2; omega += 0.011)
    {
      x.push_back(omega);
      y.push_back(proc.ReT1(1, 1, omega, k, 0, 0));
    }

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

  tLgTOtRH_massless_gluon proc(T, T * T, s1, s2, s3, s4, m1, m2, m3, m4);

  for (double k : {0.5, 1.0, 1.5})
  {
    std::vector<double> x, y;

    for (double omega = 0; omega < 2; omega += 0.0011)
    {
      x.push_back(omega);
      y.push_back(proc.ImT1(1, 1, omega, k, 0, 0));
    }

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

  tLgTOtRH_massless_gluon proc(T, T * T, s1, s2, s3, s4, m1, m2, m3, m4);

  for (double k : {0.5, 1.0, 1.5})
  {
    std::vector<double> x, y;

    for (double omega = 0; omega < 2; omega += 0.011)
    {
      x.push_back(omega);
      y.push_back(proc.ReT2(1, 1, omega, k, 0, 0));
    }

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

  tLgTOtRH_massless_gluon proc(T, T * T, s1, s2, s3, s4, m1, m2, m3, m4);

  for (double k : {0.5, 1.0, 1.5})
  {
    std::vector<double> x, y;

    for (double omega = 0; omega < 2; omega += 0.001)
    {
      x.push_back(omega);
      y.push_back(proc.ImT2(1, 1, omega, k, 0, 0));
    }

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
  test_ReT1();
  test_ImT1();
  test_ReT2();
  test_ImT2();
}

int main()
{

  // testing();

  /*int ncores = 5, pcores = 1e3;
  cubacores(&ncores, &pcores);*/

  double T = 10.;
  double v = 10;
  double m1, m2, m3, m4;

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

    t_L     h
     \     /
      \   /
       \ /
        |
        |
        |
       & \
      &   \
     &     \
    g       t_R
  */

  // Complete process. No approximations
  // Time : 8.824 s
  // Time redux :
  // T = 10 | Gamma_y = 2.92322e-06 +- 2.66048e-08
  // T = 100 | Gamma_y = 0.0132593 +- 0.000115383
  // T = 200 | Gamma_y = 0.00194427 +- 1.01971e-05
  m1 = mt; // mt
  m2 = mg; // mg
  m3 = ms; // ms
  m4 = mt; // mt
  // tLgTOtRH proc(T, T * T, s1, s2, s3, s4, m1, m2, m3, m4);

  // Complete process. No approximations but massless gluons
  // T = 10 | Gamma_y = 9.34226e-07 +- 1.26619e-08
  // T = 100 | Gamma_y = 0.000933403 +- 8.98317e-06
  // T = 200 | Gamma_y = 0.00185429 +- 1.41573e-05
  // Time : 8.037 s

  m1 = mt; // mt
  m2 = 0;  // mg
  m3 = ms; // ms
  m4 = mt; // mt
  // tLgTOtRH_massless_gluon proc(T, T * T, s1, s2, s3, s4, m1, m2, m3, m4);

  // t-channel with massless approximationsi. // not a care about helicities
  // Gamma_y = 0.00380921 +- 2.20009e-05
  // Time : 35.298 s
  // Gamma_y = 0.0376673 +- 0.000164028
  // Time redux : 34.927 s
  m1 = 0; // mt
  m2 = 0; // mg
  m3 = 0; // ms
  m4 = 0; // mt
  // tLgTOtRH_massless proc(T, T * T, s1, s2, s3, s4, m1, m2, m3, m4);

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

  /*double mA              = 10;
  double mB              = 30;
  std::vector<double> pA = {0, 0, 10};

  double phi    = 1.6;
  float cos_phi = cos(phi);
  float sin_phi = sin(phi);

  std::vector<double> pB = {10 * cos_phi, 10 * sin_phi, 0};

  std::cout << "Energy A\t" << proc.Energy(mA, pA) << "\t"
            << proc.Energy(mB, pB) << "\t" << pA * pB << "\n";

  std::cout << proc.MonteCarloInt(
                   proc.Energy(mA, pA), proc.Energy(mB, pB), pA, pB)
            << "\t\t<<<<<<<\n\n";

  // exit(0);*/
  //  t-channel with massless HTL. // care very much about helicities
  //  Gamma_y =
  //  Time :
  //  Gamma_y =
  //  Time redux :
  m1 = 0; // mt
  m2 = 0; // mg
  m3 = 0; // ms
  m4 = 0; // mt
  tLgTOtRH_massless_helicity_HTL proc(T, T * T, s1, s2, s3, s4, m1, m2, m3, m4);

  /*std::vector<double> p1 = {
      -2.4319545096067867185, 0.81642954775194209738, -2.168603517231575406};
  std::vector<double> p2 = {
      -2.9883572349034142057, 0.7114617547144796994, 4.9419805021096898656};
  double theta = M_PI;

  proc.E1_   = proc.Energy(m1, p1);
  proc.ET_   = proc.Energy(m1, p1) + proc.Energy(m2, p2);
  proc.p1_   = p1;
  proc.p2_   = p2;
  proc.p1p2_ = p1 + p2;
  proc.p12_  = proc.p1p2_ * proc.p1p2_;

  proc.integrate_phi(theta);
  exit(0);*/
  /****** useless ******/

  // identity proc(100, 1, 0, 0, 0, 0, 0, 0, 0, 0);
  // identity proc(100, 1, s1, s2, s3, s4, 0, 0, 0, 0);

  int comp, nregions, neval, fail;
  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

  auto start = std::chrono::system_clock::now();
  int mode   = 1;
  if (mode == 0)
  {
    warm_up_vegas(proc.Integrand, proc, 1e3, 20, -1, integral, error, prob);
    gridded_vegas(proc.Integrand, proc, 1e4, 10, 1, integral, error, prob);
  }
  else
  {
    Vegas(NDIM,
          NCOMP,
          proc.Integrand,
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
