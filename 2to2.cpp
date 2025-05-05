#include "2to2collisionwidth.hpp"
const int NDIM        = 6;
const int NCOMP       = 1;
const int NVEC        = 1;
const double EPSREL   = 1e-3;
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
        (pow(el, 2) * pow(gs, 2) * pow(mt, 2) *
         (2 * pow(mg, 4) * pow(ms, 2) * (pow(s, 2) + pow(t, 2)) -
          2 * pow(mg, 2) * s * t *
              (2 * pow(ms, 4) - 2 * pow(ms, 2) * (s + t) + pow(s + t, 2)) +
          pow(mt, 2) *
              (-(s * t * (pow(s, 2) + pow(t, 2))) -
               2 * pow(mg, 4) *
                   (pow(s, 2) + pow(t, 2) + 2 * pow(ms, 2) * (s + t)) +
               pow(mg, 2) * (4 * pow(ms, 4) * (s + t) + 8 * s * t * (s + t) -
                             5 * pow(ms, 2) * (pow(s, 2) + pow(t, 2))) +
               pow(mt, 2) *
                   (pow(s + t, 3) + 4 * pow(mg, 4) * (pow(ms, 2) + s + t) +
                    pow(mt, 2) *
                        (-4 * pow(mg, 4) + 6 * pow(mg, 2) * pow(ms, 2) -
                         2 * pow(mt, 4) + 4 * pow(mt, 2) * (s + t) -
                         3 * pow(s + t, 2)) -
                    2 * pow(mg, 2) *
                        (2 * pow(ms, 4) + pow(ms, 2) * (s + t) -
                         2 * (pow(s, 2) - 4 * s * t + pow(t, 2))))))) /
        (pow(mg, 2) * pow(mW, 2) * pow(-pow(mt, 2) + s, 2) * pow(sW, 2) *
         pow(-pow(mt, 2) + t, 2)));
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
        (2 * pow(el, 2) * pow(gs, 2) * pow(mt, 2) *
         (2 * pow(ms, 2) * pow(mt, 6) +
          pow(mt, 4) * (-2 * pow(ms, 4) + (s - 3 * t) * (3 * s - t)) -
          s * t * (2 * pow(ms, 4) - 2 * pow(ms, 2) * (s + t) + pow(s + t, 2)) +
          pow(mt, 2) * (2 * pow(ms, 4) * (s + t) + 4 * s * t * (s + t) -
                        3 * pow(ms, 2) * (pow(s, 2) + pow(t, 2))))) /
        (pow(mW, 2) * pow(-pow(mt, 2) + s, 2) * pow(sW, 2) *
         pow(-pow(mt, 2) + t, 2)));
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

struct tLgTOtRH_massless_helicity : Process
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

int main()
{
  /*int ncores = 5, pcores = 1e3;
  cubacores(&ncores, &pcores);*/

  double T = 100.;

  double m1, m2, m3, m4;

  double mt = 170;
  double mg = 10;
  double ms = 100;

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

    t_L     j
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
  m1 = mt; // mt
  m2 = mg; // mg
  m3 = ms; // ms
  m4 = mt; // mt
  // tLgTOtRH proc(T, T * T, s1, s2, s3, s4, m1, m2, m3, m4);

  // Complete process. No approximations but massless gluons
  // Time : 8.037 s

  m1 = mt; // mt
  m2 = 0;  // mg
  m3 = ms; // ms
  m4 = mt; // mt
  // tLgTOtRH_massless_gluon proc(T, T * T, s1, s2, s3, s4, m1, 0, m3, m4);

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
  tLgTOtRH_massless_helicity proc(T, T * T, s1, s2, s3, s4, m1, m2, m3, m4);

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
  std::cout << "\n\nGamma_y = \t" << integral[0] / (N1 * proc.prefactor)
            << " +- " << error[0] / (N1 * proc.prefactor) << "\n\n";

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
