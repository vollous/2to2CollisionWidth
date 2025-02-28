#include "2to2collisionwidth.hpp"

using namespace TwoToTwoCollisionWidth;

double Process::integrate_phi(const double &theta) {
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

  auto ptr = [=](const double &phi) -> double {
    return integrand_theta_phi(theta, phi);
  };
  gsl_function_pp<decltype(ptr)> Fp(ptr);
  const gsl_function *F = static_cast<gsl_function *>(&Fp);

  double result, error;
  gsl_integration_qag(F, 0, M_PI, 1e-7, 1e-7, 1000, GSL_INTEG_GAUSS15, w,
                      &result, &error);
  gsl_integration_workspace_free(w);
  return result;
}

double Process::Integrate(const double &E1, const double &ET, const double &s,
                          const std::vector<double> &p1,
                          const std::vector<double> &p1p2) {
  E1_ = E1;
  ET_ = ET;
  s_ = s;
  p1_ = p1;
  p1p2_ = p1p2;
  p12_ = p1p2_ * p1p2_;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
  gsl_function F;

  // Capture 'this' using another lambda
  auto lambda = [](double theta, void *params) -> double {
    return static_cast<Process *>(params)->integrate_phi(theta);
  };

  F.function = lambda;
  F.params = this;

  double result, error;
  gsl_integration_qag(&F, 0, 2 * M_PI, 1e-7, 1e-7, 1000, GSL_INTEG_GAUSS15, w,
                      &result, &error);
  gsl_integration_workspace_free(w);
  return result;
}

double Process::Energy(const double &m, const std::vector<double> &p) {
  return sqrt(pow(m, 2) + pow(p[0], 2) + pow(p[1], 2) + pow(p[2], 2));
}

double Process::Distribution(const double &m, const std::vector<double> &p,
                             const int &s) {
  return 1 / (exp(Energy(m, p) / T) + s);
}

void Process::calculate_r(double &r, double &delta_r, double shift,
                          double theta, double phi, const double &rcoeff) {

  double sqrtpart = sqrt(
      pow(rcoeff, 2) * pow(-4 * pow(ET_, 2) + 4 * pow(m3, 2) - 4 * pow(m4, 2) +
                               4 * p1p2_ * p1p2_ + 8 * pow(ET_, 2) * shift -
                               8 * p1p2_ * p1p2_ * shift,
                           2) -
      4 * (4 * pow(ET_, 2) - 4 * pow(rcoeff, 2)) *
          (-pow(ET_, 4) + 2 * pow(ET_, 2) * pow(m3, 2) - pow(m3, 4) +
           2 * pow(ET_, 2) * pow(m4, 2) + 2 * pow(m3, 2) * pow(m4, 2) -
           pow(m4, 4) + 2 * pow(ET_, 2) * p1p2_ * p1p2_ -
           2 * pow(m3, 2) * p1p2_ * p1p2_ + 2 * pow(m4, 2) * p1p2_ * p1p2_ -
           pow(p1p2_ * p1p2_, 2) - 4 * pow(ET_, 2) * p1p2_ * p1p2_ * shift +
           4 * pow(m3, 2) * p1p2_ * p1p2_ * shift -
           4 * pow(m4, 2) * p1p2_ * p1p2_ * shift +
           4 * pow(p1p2_ * p1p2_, 2) * shift +
           4 * pow(ET_, 2) * p1p2_ * p1p2_ * pow(shift, 2) -
           4 * pow(p1p2_ * p1p2_, 2) * pow(shift, 2)));

  if (sqrtpart < 0) {
    std::cout << "Check\t"
              << (pow(4 * pow(m3, 2) - 4 * pow(m4, 2), 2) * pow(rcoeff, 2) -
                  4 *
                      (-4 * pow(ET_, 4) + 8 * pow(ET_, 2) * pow(m3, 2) -
                       4 * pow(m3, 4) + 8 * pow(ET_, 2) * pow(m4, 2) +
                       8 * pow(m3, 2) * pow(m4, 2) - 4 * pow(m4, 4) +
                       4 * pow(ET_, 2) * p12_) *
                      (16 * pow(ET_, 2) - pow(rcoeff, 2)))
              << "\n";
    std::cout << "\nCannot calculate rr\n\n";
    std::cout << "shift\t" << shift << "\n";
    std::cout << "p1_\t" << p1_[0] << "\t" << p1_[1] << "\t" << p1_[2] << "\t"
              << "\n";
    std::cout << "p2_\t" << (p1p2_ - p1_)[0] << "\t" << (p1p2_ - p1_)[1] << "\t"
              << (p1p2_ - p1_)[2] << "\t"
              << "\n";
    std::cout << "p1 + p2\t" << p1p2_[0] << "\t" << p1p2_[1] << "\t" << p1p2_[2]
              << "\n";
    std::cout << "p1 . p2\t" << p12_ << "\n";
    std::cout << "theta\t" << theta << "\n";
    std::cout << "phi\t" << phi << "\n";
    std::cout << "rcoeff\t" << rcoeff << "\n";
    std::cout << "Spare energy\t" << ET_ - Energy(m3 + m4, p1p2_) << "\n";
    std::cout << "p12_\t" << p12_ << "\n";
    std::cout << "ET_\t" << ET_ << "\n";
    std::cout << "m3\t" << m3 << "\n";
    std::cout << "m4\t" << m4 << "\n";

    throw;
  };

  r = (-(rcoeff * (-4 * pow(ET_, 2) + 4 * pow(m3, 2) - 4 * pow(m4, 2) +
                   4 * p1p2_ * p1p2_ + 8 * pow(ET_, 2) * shift -
                   8 * p1p2_ * p1p2_ * shift)) +
       sqrtpart) /
      (2. * (4 * pow(ET_, 2) - 4 * pow(rcoeff, 2)));

  delta_r =
      abs(-0.5 * (2 * r + 2 * rcoeff * (-1 + shift)) /
              sqrt(pow(m3, 2) + pow(r, 2) + 2 * rcoeff * r * (-1 + shift) +
                   p1p2_ * p1p2_ * pow(-1 + shift, 2)) -
          (2 * r + 2 * rcoeff * shift) /
              (2. * sqrt(pow(m4, 2) + pow(r, 2) + 2 * rcoeff * r * shift +
                         p1p2_ * p1p2_ * pow(shift, 2))));
}

double Process::integrand_theta_phi(const double &theta, const double &phi) {

  const std::vector<double> p4_centered = {cos(theta) * sin(phi),
                                           sin(theta) * sin(phi), cos(phi)};

  const double rcoeff = p4_centered * p1p2_;
  double r, delta_r;
  const double shift = (pow(ET_, 2) - pow(m3, 2) + pow(m4, 2) - p1p2_ * p1p2_) /
                       (2 * pow(ET_, 2) - 2 * p1p2_ * p1p2_);
  calculate_r(r, delta_r, shift, theta, phi, rcoeff);
  if (r < 0) {
    std::cout << "Negative r\n";
    throw;
  }

  const std::vector<double> p4 = shift * p1p2_ + r * p4_centered;
  const std::vector<double> p3 = p1p2_ - p4;

  if (abs(ET_ / (Energy(m3, p3) + Energy(m4, p4))) - 1 > 1e-8) {
    std::cout << "------------------\tEnergy mismatch\t---------------------\n";
    std::cout << "r\t" << r << "\n";
    std::cout << "shift\t" << shift << "\n";
    std::cout << "rcoeff\t" << rcoeff << "\n";
    std::cout << "p12\t" << p12_ << "\n";
    std::cout << "delta_r\t" << delta_r << "\n";
    std::cout << "p1\t" << p1_[0] << "\t" << p1_[1] << "\t" << p1_[2] << "\n";
    std::cout << "p2\t" << (p1p2_ - p1_)[0] << "\t" << (p1p2_ - p1_)[1] << "\t"
              << (p1p2_ - p1_)[2] << "\n";
    std::cout << "p3\t" << p3[0] << "\t" << p3[1] << "\t" << p3[2] << "\n";
    std::cout << "p4\t" << p4[0] << "\t" << p4[1] << "\t" << p4[2] << "\n";

    std::cout << "\n--------\n";

    std::cout << "E1\t" << E1_ << "\n";
    std::cout << "E2\t" << ET_ - E1_ << "\n";
    std::cout << "E3\t" << Energy(m3, p3) << "\n";
    std::cout << "E4\t" << Energy(m4, p4) << "\n";

    std::cout << "\nE1 + E2\t" << ET_ << "\n";
    std::cout << "E3 + E4\t" << Energy(m3, p3) + Energy(m4, p4) << "\n";

    std::cout << "--------\n";

    std::cout << "m1\t" << m1 << "\n";
    std::cout << "m2\t" << m2 << "\n";
    std::cout << "m3\t" << m3 << "\n";
    std::cout << "m4\t" << m4 << "\n";

    std::cout << "--------\n";

    std::cout << "theta\t" << theta << "\n";
    std::cout << "phi\t" << phi << "\n";

    throw;
  }

  const double t =
      pow(m1, 2) + pow(m3, 2) - 2 * E1_ * Energy(m3, p3) + p1_ * p3;

  return AmplitudeSquared(s_, t) * (1 - Distribution(m3, p3, s3)) *
         (1 + Distribution(m4, p4, s4)) * r * r * sin(phi) /
         (delta_r * 4 * Energy(m3, p3) * Energy(m4, p4));
}

double Process::MonteCarloInt(const std::vector<double> p1,
                              const std::vector<double> p2) {

  const std::vector<double> p1p2 = p1 + p2;
  const double p1dotp2 = p1 * p2;
  const double E1 = Energy(m1, p1);
  const double E2 = Energy(m2, p2);
  const double s = m1 * m1 + m2 * m2 + 2 * E1 * E2 - 2 * p1dotp2;

  // Theta function to see if its possible to have a
  // process with a given p1 and p2
  if (E1 + E2 - Energy(m3 + m4, p1p2) <= 0)
    return 0.;

  double integral = Integrate(E1, E1 + E2, s, p1, p1p2);

  // Normalization
  integral *=
      Distribution(m1, p1, s1) * Distribution(m2, p2, s2) / (4 * E1 * E2);

  return integral;
};

Process::Process(const double &T_in, const int &s1_in, const int &s2_in,
                 const int &s3_in, const int &s4_in, const double &m1_in,
                 const double &m2_in, const double &m3_in, const double &m4_in)
    : s1(s1_in), s2(s2_in), s3(s3_in), s4(s4_in), T(T_in), m1(m1_in), m2(m2_in),
      m3(m3_in), m4(m4_in) {};

struct tLgTOtRH : Process {
  using Process::Process; // Import constructor

  double AmplitudeSquared(const double &s, const double &t) override {
    const double e = 0.31;
    const double gs = 1.;
    const double sw = 0.472;
    const double mw = 80;

    return (2 * pow(e, 2) * pow(gs, 2) * pow(m1, 2) *
            (-(pow(m1, 6) *
               (pow(m2, 4) + pow(m2, 2) * (-3 * pow(m3, 2) +
                                           17 * (2 * pow(m1, 2) + pow(m2, 2) +
                                                 pow(m3, 2) - s - t)))) +
             pow(m1, 4) *
                 (pow(m2, 4) *
                      (2 * pow(m3, 2) + s -
                       4 * (2 * pow(m1, 2) + pow(m2, 2) + pow(m3, 2) - s - t) +
                       4 * t) +
                  pow(m2, 2) *

                      (-3 * pow(m3, 4) +
                       pow(m3, 2) * (-9 * pow(m3, 2) +
                                     4 * (2 * pow(m1, 2) + pow(m2, 2) +
                                          pow(m3, 2) - s - t) -
                                     6 * t) +
                       pow(m3, 2) * (10 * (2 * pow(m1, 2) + pow(m2, 2) +
                                           pow(m3, 2) - s - t) +
                                     t) +
                       (2 * pow(m1, 2) + pow(m2, 2) + pow(m3, 2) - s - t) *
                           (2 * (2 * pow(m1, 2) + pow(m2, 2) + pow(m3, 2) - s -
                                 t) +
                            31 * t))) +
             t * (-(pow(m2, 4) * s * (2 * pow(m3, 2) - 3 * t)) -
                  pow(m2, 2) *
                      (pow(m3, 6) + pow(m3, 4) * t -
                       pow(m3, 2) *
                           (2 * t *
                                (3 * (2 * pow(m1, 2) + pow(m2, 2) + pow(m3, 2) -
                                      s - t) +
                                 t) +
                            pow(m3, 2) * (-2 * pow(m1, 2) - pow(m2, 2) -
                                          pow(m3, 2) + s + 2 * t)) +
                       (2 * pow(m1, 2) + pow(m2, 2) + pow(m3, 2) - s) *
                           (pow(m3, 2) * (-2 * pow(m1, 2) - pow(m2, 2) -
                                          pow(m3, 2) + s + 2 * t) +
                            (2 * pow(m1, 2) + pow(m2, 2) + pow(m3, 2) - s - t) *
                                (-2 * pow(m1, 2) - pow(m2, 2) - pow(m3, 2) + s +
                                 4 * t)))) +
             pow(m1, 2) *
                 (pow(m2, 4) *
                      (4 * s *
                           (2 * pow(m1, 2) + pow(m2, 2) + pow(m3, 2) - s - t) -
                       4 * s * t - 3 * pow(t, 2) + pow(m3, 2) * (-3 * s + t) +
                       pow(m3, 2) * (s + t)) +
                  pow(m2, 2) *
                      (5 * pow(m3, 6) -
                       4 * pow(m3, 4) *
                           (2 * (2 * pow(m1, 2) + pow(m2, 2) + pow(m3, 2) - s -
                                 t) +
                            3 * t) +
                       pow(m3, 2) * (-pow(2 * pow(m1, 2) + pow(m2, 2) +
                                              pow(m3, 2) - s - t,
                                          2) -
                                     10 *
                                         (2 * pow(m1, 2) + pow(m2, 2) +
                                          pow(m3, 2) - s - t) *
                                         t +
                                     pow(t, 2)) -
                       (2 * pow(m1, 2) + pow(m2, 2) + pow(m3, 2) - s - t) *
                           (pow(2 * pow(m1, 2) + pow(m2, 2) + pow(m3, 2) - s -
                                    t,
                                2) +
                            16 *
                                (2 * pow(m1, 2) + pow(m2, 2) + pow(m3, 2) - s -
                                 t) *
                                t +
                            11 * pow(t, 2)) +
                       pow(m3, 2) *
                           (pow(m3, 2) * (2 * pow(m1, 2) + pow(m2, 2) +
                                          pow(m3, 2) - s + 7 * t) +
                            2 * (2 * pow(m1, 2) + pow(m2, 2) + pow(m3, 2) - s - t) *
                                (2 * (2 * pow(m1, 2) + pow(m2, 2) + pow(m3, 2) -
                                      s - t) +
                                 11 * t)))))) /
           (pow(m2, 2) * pow(mw, 2) * pow(-pow(m1, 2) + s, 2) * pow(sw, 2) *
            pow(-pow(m1, 2) + t, 2));
  }
};

static int Integrand(const int *ndim, const cubareal xx[], const int *ncomp,
                     cubareal ff[], void *userdata) {
  auto proc = static_cast<tLgTOtRH *>(userdata);

  if (xx[0] == 1 or xx[3] == 1) {
    // TODO might be removable
    ff[0] = 0.;
    return 0;
  }

  const double scalling = 100; // Put the hypercube on the region of interest.

  // Spherical coordinates p1
  const double r1 = scalling * xx[0] / (1 - xx[0]);
  const double theta1 = 2 * M_PI * xx[1];
  const double phi1 = M_PI * xx[2];
  const std::vector<double> p1 =
      r1 * (std::vector<double>){cos(theta1) * sin(phi1),
                                 sin(theta1) * sin(phi1), cos(phi1)};

  // Spherical coordinates p2
  const double r2 = scalling * xx[3] / (1 - xx[3]);
  const double theta2 = 2 * M_PI * xx[4];
  const double phi2 = M_PI * xx[5];
  const std::vector<double> p2 =
      r2 * (std::vector<double>){cos(theta2) * sin(phi2),
                                 sin(theta2) * sin(phi2), cos(phi2)};
  // Return result
  // std::cout << "p1 >\t" << p1[0] << "\t" << p1[1] << "\t" << p1[2] << "\t";
  // std::cout << "p2 >\t" << p2[0] << "\t" << p2[1] << "\t" << p2[2] << "\t";

  ff[0] = proc->MonteCarloInt(p1, p2) * _4_M_4 * pow(r1, 2) *
          pow(r1 + scalling, 2) * sin(phi1) * pow(r2, 2) * sin(phi2) *
          pow(r2 + scalling, 2) / (scalling * scalling);

  // std::cout << "r >\t" << ff[0] << "\n";
  return 0;
};

int main() {
  double T = 100.;

  int s1 = 1;
  int s2 = -1;
  int s3 = 1;
  int s4 = -1;

  double m1 = 170;
  double m2 = 10;
  double m3 = 170;
  double m4 = 100;

  tLgTOtRH proc(T, s1, s2, s3, s4, m1, m2, m3, m4);

  const int NDIM = 6;
  const int NCOMP = 1;
  const int NVEC = 1;
  const double EPSREL = 5e-2;
  const double EPSABS = 1e-12;
  const int VERBOSE = 3;
  const int LAST = 0;
  const int SEED = 0;
  const int MINEVAL = 0;
  const int MAXEVAL = 200000;
  const int NSTART = 1000;
  const int NINCREASE = 1000;
  const int NBATCH = 1000;
  const int GRIDNO = 0;
  const char *STATEFILE = NULL;
  const int *SPIN = NULL;

  int comp, nregions, neval, fail;
  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

  Vegas(NDIM, NCOMP, Integrand, &proc, NVEC, EPSREL, EPSABS, VERBOSE, SEED,
        MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH, GRIDNO, STATEFILE, &SPIN,
        &neval, &fail, integral, error, prob);

  std::cout << "Integral\t" << integral[0] << "\n";
  std::cout << "Error \t" << error[0] << "\n";
  std::cout << "Relative error \t" << error[0] / integral[0] << "\n";
  std::cout << "prob\t" << prob[0] << "\n";

  std::cout << "\n\nGamma_y = \t" << integral[0] * 12 / pow(T, 3) << "\n";
  return 0;
}