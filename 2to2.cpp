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
using namespace TwoToTwoCollisionWidth;

struct tLgTOtRH : Process
{
  using Process::Process; // Import constructor

  double AmplitudeSquared(const double &s, const double &t) override
  {
    const double e  = 0.31;
    const double gs = 1.;
    const double sw = 0.472;
    const double mw = 80;

    return (pow(e, 2) * pow(gs, 2) * pow(m1, 4) *
            pow(-2 * pow(m1, 2) - pow(m2, 2) + s + t, 2) *
            (-pow(m1, 4) - pow(m2, 4) +
             pow(m2, 2) * (pow(m3, 2) - 2 * (2 * pow(m1, 2) + pow(m2, 2) +
                                             pow(m3, 2) - s - t)) -
             s * t + pow(m1, 2) * (2 * pow(m2, 2) + s + t))) /
           (pow(m2, 2) * pow(mw, 2) * pow(-pow(m1, 2) + s, 2) * pow(sw, 2) *
            pow(-pow(m1, 2) + t, 2));
  }
};

struct identity : Process
{
  using Process::Process; // Import constructor

  double AmplitudeSquared(const double &s, const double &t) override
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
  // identity proc(100, 0, 0, 0, 0, 0, 0, 0, 0);
  identity proc(100, s1, s2, s3, s4, 0, 0, 0, 0);

  int comp, nregions, neval, fail;
  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

  auto start = std::chrono::system_clock::now();
  int mode   = 0;
  if (mode == 0)
  {
    warm_up_vegas(proc.Integrand, proc, 1e4, 20, -1, integral, error, prob);
    gridded_vegas(proc.Integrand, proc, 1e5, 10, 1, integral, error, prob);
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

  std::cout << "\n\nGamma_y = \t" << integral[0] * 12. / pow(T, 3) << " +- "
            << error[0] * 12. / pow(T, 3) << "\n\n";

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
