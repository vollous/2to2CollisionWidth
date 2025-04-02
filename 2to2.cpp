#include "2to2collisionwidth.hpp"
const int NDIM        = 6;
const int NCOMP       = 1;
const int NVEC        = 1;
const double EPSREL   = 5e-2;
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
    const double s =
        m1 * m1 + m2 * m2 + 2 * Energy(m1, p1) * Energy(m2, p2) - 2 * p1dotp2;

    return (pow(el, 2) * pow(gs, 2) * pow(m1, 4) *
            pow(-2 * pow(m1, 2) - pow(m2, 2) + s + t, 2) *
            (-pow(m1, 4) - pow(m2, 4) +
             pow(m2, 2) * (pow(m3, 2) - 2 * (2 * pow(m1, 2) + pow(m2, 2) +
                                             pow(m3, 2) - s - t)) -
             s * t + pow(m1, 2) * (2 * pow(m2, 2) + s + t))) /
           (pow(m2, 2) * pow(mW, 2) * pow(-pow(m1, 2) + s, 2) * pow(sW, 2) *
            pow(-pow(m1, 2) + t, 2));
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
  double T = 100.;

  int s1 = 1;
  int s2 = -1;
  int s3 = 1;
  int s4 = -1;

  double m1 = 170;
  double m2 = 10;
  double m3 = 170;
  double m4 = 100;

  // tLgTOtRH proc(T, s1, s2, s3, s4, m1, m2, m3, m4);
  tLgTOtRH_massless proc(T, s1, s2, s3, s4, 0, 0, 0, 0);
  // identity proc(100, 0, 0, 0, 0, 0, 0, 0, 0);
  // identity proc(100, s1, s2, s3, s4, 0, 0, 0, 0);

  int comp, nregions, neval, fail;
  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

  auto start = std::chrono::system_clock::now();
  int mode   = 0;
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

  const double vw            = 0.95;
  const double gamma         = 1; /*/ sqrt(1 - pow(vw, 2));*/
  const double N1            = gamma * 2 * pow(M_PI, 3) * pow(T, 2) / 3;
  const double TempPrefactor = T * T;
  std::cout << "\n\nGamma_y = \t" << integral[0] / (N1 * TempPrefactor)
            << " +- " << error[0] / (N1 * TempPrefactor) << "\n\n";

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
