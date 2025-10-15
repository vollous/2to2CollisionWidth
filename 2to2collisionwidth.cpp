#include "2to2collisionwidth.hpp"

using namespace TwoToTwoCollisionWidth;

void my_gsl_error_handler(const char *reason,
                          const char *file,
                          int line,
                          int gsl_errno)
{
  std::cerr << "GSL Error: " << reason << " in " << file << ":" << line
            << " (Error Code: " << gsl_errno << ")" << std::endl;

  return;
}

double Process::integrate_phi(const double &theta)
{
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);

  gsl_set_error_handler(&my_gsl_error_handler);

  auto ptr = [=](const double &phi) -> double
  { return integrand_theta_phi(theta, phi); };
  gsl_function_pp<decltype(ptr)> Fp(ptr);
  const gsl_function *F = static_cast<gsl_function *>(&Fp);

  double result, error;
  if (gsl_integration_qag(F,
                          0,
                          M_PI,
                          1e-2,
                          1e-2,
                          10000,
                          GSL_INTEG_GAUSS61,
                          w,
                          &result,
                          &error) != 0)
  {
    std::cout << "error\t" << error << "\n";
    std::cout << "theta\t" << theta << "\n";
    std::cout << "theta - pi\t" << std::setprecision(20) << theta - M_PI
              << "\n";

    for (double p = 0; p < M_PI; p += 0.01)
    {
      std::cout << p << "\t" << ptr(p) << "\n";
    }

    std::cout << "E1\t" << E1_ << "\n";
    std::cout << "E2\t" << ET_ - E1_ << "\n";

    std::cout << "\nE1 + E2\t" << ET_ << "\n";
    std::cout << "Minimum Energy\t" << (ET_ - Energy(m3 + m4, p1p2_)) << "\n";
    exit(0);
  }
  gsl_integration_workspace_free(w);
  return result;
}

double Process::Integrate(const double &E1,
                          const double &ET,
                          const std::vector<double> &p1,
                          const std::vector<double> &p2,
                          const std::vector<double> &p1p2)
{
  E1_   = E1;
  ET_   = ET;
  p1_   = p1;
  p2_   = p2;
  p1p2_ = p1p2;
  p12_  = p1p2_ * p1p2_;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);
  gsl_set_error_handler(&my_gsl_error_handler);
  gsl_function F;

  // Capture 'this' using another lambda
  auto lambda = [](double theta, void *params) -> double
  { return static_cast<Process *>(params)->integrate_phi(theta); };

  F.function = lambda;
  F.params   = this;

  double result, error;

  gsl_integration_qag(&F,
                      0,
                      2 * M_PI,
                      1e-2,
                      1e-2,
                      10000,
                      GSL_INTEG_GAUSS61,
                      w,
                      &result,
                      &error);
  gsl_integration_workspace_free(w);
  return result;
}

double Process::Energy(const double &m, const std::vector<double> &p)
{
  return sqrt(pow(m, 2) + pow(p[0], 2) + pow(p[1], 2) + pow(p[2], 2));
}

double Process::Distribution(const double &m,
                             const std::vector<double> &p,
                             const int &s)
{
  return 1 / (exp(Energy(m, p) / T) + s);
}

double Process::HTLa(const double &m, const double &omega, const double &k)
{
  const double logomegak = log(abs((omega + k) / (omega - k)));
  return pow(m, 2) / k * (1 - omega / (2. * k) * logomegak);
}

double Process::HTLb(const double &m, const double &omega, const double &k)
{
  const double logomegak = log(abs((omega + k) / (omega - k)));
  return pow(m, 2) / k *
         (-omega / k + (pow(omega / k, 2) - 1) / 2. * logomegak);
}

double Process::L1(const double &p, const double &omega, const double &k)
{
  return log(abs((pow(k - 2 * p, 2) - pow(omega, 2)) /
                 (pow(k + 2 * p, 2) - pow(omega, 2))));
}

double Process::L2(const double &p, const double &omega, const double &k)
{
  return log(abs((pow(omega + k, 2) * (pow(omega - k, 2) - 4 * pow(p, 2))) /
                 (pow(omega - k, 2) * (pow(omega + k, 2) - 4 * pow(p, 2)))));
}

double Process::ReT1(const double &g,
                     const double &C,
                     const double &omega,
                     const double &k,
                     const double &mB,
                     const double &mF)
{
  auto ptr = [=](const double &p) -> double
  {
    if (p == 0 and mB == 0)
      return (-8 * k * T * omega) / (pow(k, 2) - pow(omega, 2));
    const double _L1 = L1(p, omega, k);
    const double _L2 = L2(p, omega, k);
    const double _nF = Distribution(mF, {p, 0, 0}, 1);
    const double _nB = Distribution(mB, {p, 0, 0}, -1);

    return (p * _L2 * (_nB + _nF) + omega * _L1 * _nB);
  };

  double r;

  gsl_function_pp<decltype(ptr)> Fp(ptr);
  gsl_function *F = static_cast<gsl_function *>(&Fp);

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);
  gsl_set_error_handler(&my_gsl_error_handler);

  double result, error;
  if (gsl_integration_qagiu(
          F, omega + k, 1e-6, 1e-6, 10000, w, &result, &error) != 0)
  {
    std::cout << "Real part of T1 (1))\n";
    std::cout << "error\t" << error << "\n";
    std::cout << "error\t" << error << "\n";
    std::cout << "omega\t" << omega << "\n";
    std::cout << "k\t" << k << "\n";
    std::cout << "omega - k\t" << omega - k << "\n";
    exit(0);
  }
  gsl_integration_workspace_free(w);
  r += result;

  gsl_integration_workspace *w1 = gsl_integration_workspace_alloc(10000);
  size_t np;
  double *points_ptr;
  double points_zero[3]{0., (omega + k) / 2., omega + k};
  double points[4]{0., abs(omega - k) / 2., (omega + k) / 2., omega + k};
  if (omega == 0 or k == 0)
  {
    np         = 3;
    points_ptr = points_zero;
  }
  else
  {
    np         = 4;
    points_ptr = points;
  }

  double result1, error1;
  if (gsl_integration_qagp(
          F, points_ptr, np, 1e-6, 1e-6, 10000, w1, &result1, &error1) != 0)
  {
    std::cout << "Real part of T1 (2))\n";
    std::cout << "result1\t" << result1 << "\n";
    std::cout << "error1\t" << error1 << "\n";
    std::cout << "omega\t" << omega << "\n";
    std::cout << "k\t" << k << "\n";
    std::cout << "omega - k\t" << omega - k << "\n";
    exit(0);
    // return 0;
  }
  gsl_integration_workspace_free(w1);

  r += result1;

  return g * g * C * r / (8. * M_PI * M_PI * k);
}
double Process::ImT1(const double &g,
                     const double &C,
                     const double &omega,
                     const double &k,
                     const double &mB,
                     const double &mF)
{
  if (abs(omega - k) < 1e-10) return -100.;

  double r = 0;

  if (omega > k)
  {
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);

    gsl_set_error_handler(&my_gsl_error_handler);

    auto ptr = [=](const double &p) -> double
    {
      return (Distribution(mB, {p, 0, 0}, -1) * (p - omega) -
              p * Distribution(mF, {p, 0, 0}, 1));
    };

    gsl_function_pp<decltype(ptr)> Fp(ptr);
    const gsl_function *F = static_cast<gsl_function *>(&Fp);

    double result, error;
    if (gsl_integration_qag(F,
                            abs(omega - k) / 2.,
                            (omega + k) / 2.,
                            1e-2,
                            1e-2,
                            10000,
                            GSL_INTEG_GAUSS61,
                            w,
                            &result,
                            &error) != 0)
    {
      std::cout << "Imaginary part of T1 (1))\n";
      std::cout << "error\t" << error << "\n";
      std::cout << "omega\t" << omega << "\n";
      std::cout << "k\t" << k << "\n";
      std::cout << "omega - k\t" << omega - k << "\n";
      exit(0);
    }
    gsl_integration_workspace_free(w);
    r += result;
  }
  else
  {

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);

    gsl_set_error_handler(&my_gsl_error_handler);

    auto ptr = [=](const double &p) -> double
    {
      return -(Distribution(mB, {p, 0, 0}, -1) * (p - omega) +
               p * Distribution(mF, {p, 0, 0}, 1));
    };

    gsl_function_pp<decltype(ptr)> Fp(ptr);
    gsl_function *F = static_cast<gsl_function *>(&Fp);

    double result, error;
    if (gsl_integration_qagiu(
            F, (omega + k) / 2., 1e-6, 1e-6, 10000, w, &result, &error) != 0)
    {
      std::cout << "Imaginary part of T1 (2))\n";
      std::cout << "error\t" << error << "\n";
      std::cout << "omega\t" << omega << "\n";
      std::cout << "k\t" << k << "\n";
      std::cout << "omega - k\t" << omega - k << "\n";
      exit(0);
    }
    gsl_integration_workspace_free(w);
    gsl_integration_workspace *w2 = gsl_integration_workspace_alloc(10000);
    r += result; // This sign comes from the retarded propagator

    auto ptr2 = [=](const double &p) -> double
    {
      return -(Distribution(mB, {p, 0, 0}, -1) * (p + omega) +
               p * Distribution(mF, {p, 0, 0}, 1));
    };
    gsl_function_pp<decltype(ptr2)> Fp2(ptr2);
    gsl_function *F2 = static_cast<gsl_function *>(&Fp2);

    result, error;
    if (gsl_integration_qagiu(
            F2, (k - omega) / 2., 1e-6, 1e-6, 10000, w2, &result, &error) != 0)
    {
      std::cout << "Imaginary part of T1 (3))\n";
      std::cout << "error\t" << error << "\n";
      std::cout << "omega\t" << omega << "\n";
      std::cout << "k\t" << k << "\n";
      std::cout << "omega - k\t" << omega - k << "\n";
      exit(0);
    }
    gsl_integration_workspace_free(w);
    r += result;
  }
  return g * g * C * r / (8. * M_PI * k);
}
double Process::ReT2(const double &g,
                     const double &C,
                     const double &omega,
                     const double &k,
                     const double &mB,
                     const double &mF)
{
  const double fac = (omega * omega - k * k) / (2. * k);
  auto ptr         = [=](const double &p) -> double
  {
    if (p == 0.) return 0.;
    const double _L1 = L1(p, omega, k);
    const double _nF = Distribution(mF, {p, 0, 0}, 1);
    const double _nB = Distribution(mB, {p, 0, 0}, -1);

    return (4 * p + fac * _L1) * _nB + (4 * p - fac * _L1) * _nF;
  };

  double r;

  gsl_function_pp<decltype(ptr)> Fp(ptr);
  gsl_function *F = static_cast<gsl_function *>(&Fp);

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);
  gsl_set_error_handler(&my_gsl_error_handler);

  double result, error;
  if (gsl_integration_qagiu(
          F, omega + k, 1e-6, 1e-6, 10000, w, &result, &error) != 0)
  {
    std::cout << "Real part of T2 (1))\n";
    std::cout << "error\t" << error << "\n";
    std::cout << "error\t" << error << "\n";
    std::cout << "omega\t" << omega << "\n";
    std::cout << "k\t" << k << "\n";
    std::cout << "omega - k\t" << omega - k << "\n";
    exit(0);
  }
  gsl_integration_workspace_free(w);
  r += result;
  gsl_integration_workspace *w1 = gsl_integration_workspace_alloc(10000);

  size_t np;
  double *points_ptr;
  double points_zero[3]{0., (omega + k) / 2., omega + k};
  double points[4]{0., abs(omega - k) / 2., (omega + k) / 2., omega + k};
  if (omega == 0 or k == 0)
  {
    np         = 3;
    points_ptr = points_zero;
  }
  else
  {
    np         = 4;
    points_ptr = points;
  }

  double result1, error1;
  if (gsl_integration_qagp(
          F, points_ptr, np, 1e-6, 1e-6, 10000, w, &result1, &error1) != 0)
  {
    std::cout << "np\t" << np << "\n";
    for (int i = 0; i < sizeof(np); i++)
    {
      std::cout << ">" << i << "\t" << points_ptr[i] << "\n";
    }
    std::cout << "Real part of T2 (2))\n";
    std::cout << "result1\t" << result1 << "\n";
    std::cout << "error1\t" << error1 << "\n";
    std::cout << "omega\t" << omega << "\n";
    std::cout << "k\t" << k << "\n";
    std::cout << "omega - k\t" << omega - k << "\n";
    exit(0);
    // return 0;
  }
  gsl_integration_workspace_free(w1);

  r += result1;

  return g * g * C * r / (8. * M_PI * M_PI);
}
double Process::ImT2(const double &g,
                     const double &C,
                     const double &omega,
                     const double &k,
                     const double &mB,
                     const double &mF)
{
  if (abs(omega - k) < 1e-10) return 0.;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);

  gsl_set_error_handler(&my_gsl_error_handler);

  auto ptr = [=](const double &p) -> double
  {
    return (Distribution(mB, {p, 0, 0}, -1) - Distribution(mF, {p, 0, 0}, 1));
  };

  gsl_function_pp<decltype(ptr)> Fp(ptr);
  const gsl_function *F = static_cast<gsl_function *>(&Fp);

  double result, error;
  if (gsl_integration_qag(F,
                          abs(omega - k) / 2.,
                          (omega + k) / 2.,
                          1e-6,
                          1e-6,
                          10000,
                          GSL_INTEG_GAUSS61,
                          w,
                          &result,
                          &error) != 0)
  {
    std::cout << "Imaginary part of T2\n";
    std::cout << "error\t" << error << "\n";
    std::cout << "omega\t" << omega << "\n";
    std::cout << "k\t" << k << "\n";
    std::cout << "omega - k\t" << omega - k << "\n";
    return INFINITY;
  }
  gsl_integration_workspace_free(w);

  return -g * g * C / (16. * M_PI * k) * (omega * omega - k * k) * result;
}

void Process::calculate_r(double &r,
                          double &delta_r,
                          double shift,
                          double theta,
                          double phi,
                          const double &rcoeff)
{

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

  if (sqrtpart < 0)
  {
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

double Process::integrand_theta_phi(const double &theta, const double &phi)
{

  const std::vector<double> p4_centered = {
      cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi)};

  const double rcoeff = p4_centered * p1p2_;
  double r, delta_r;
  const double shift = (pow(ET_, 2) - pow(m3, 2) + pow(m4, 2) - p1p2_ * p1p2_) /
                       (2 * pow(ET_, 2) - 2 * p1p2_ * p1p2_);
  calculate_r(r, delta_r, shift, theta, phi, rcoeff);
  if (r < 0)
  {
    std::cout << "Negative r\n";
    throw;
  }

  const std::vector<double> p4 = shift * p1p2_ + r * p4_centered;
  const std::vector<double> p3 = p1p2_ - p4;

  if (abs(ET_ / (Energy(m3, p3) + Energy(m4, p4)) - 1) > 1e-8)
  {
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

  return AmplitudeSquared(p1_, p2_, p3) * (1 - Distribution(m3, p3, s3)) *
         (1 + Distribution(m4, p4, s4)) * r * r * sin(phi) /
         (delta_r * 4 * Energy(m3, p3) * Energy(m4, p4));
}

double Process::MonteCarloInt(const double &E1,
                              const double &E2,
                              const std::vector<double> p1,
                              const std::vector<double> p2)
{

  const std::vector<double> p1p2 = p1 + p2;

  // Theta function to see if its possible to have a
  // process with a given p1 and p2
  if (E1 + E2 - Energy(m3 + m4, p1p2) <= 0) return 0.;

  double integral = Integrate(E1, E1 + E2, p1, p2, p1p2);

  // Normalization
  integral *=
      _2_PI_FACTORS * Distribution(m1, p1, s1) * Distribution(m2, p2, s2);

  return integral;
};

Process::Process(const double &T_in,
                 const double &prefactor_in,
                 const int &s1_in,
                 const int &s2_in,
                 const int &s3_in,
                 const int &s4_in,
                 const double &m1_in,
                 const double &m2_in,
                 const double &m3_in,
                 const double &m4_in)
    : s1(s1_in)
    , s2(s2_in)
    , s3(s3_in)
    , s4(s4_in)
    , T(T_in)
    , prefactor(prefactor_in)
    , m1(m1_in)
    , m2(m2_in)
    , m3(m3_in)
    , m4(m4_in) {};

int Process::Integrand(const int *ndim,
                       const cubareal xx[],
                       const int *ncomp,
                       cubareal ff[],
                       void *userdata)
{
  auto proc = static_cast<Process *>(userdata);

  if (xx[0] == 1 or xx[3] == 1)
  {
    // TODO might be removable
    ff[0] = 0.;
    return 0;
  }

  const double scalling =
      proc->T; // Put the hypercube on the region of interest.

  // Spherical coordinates p1
  const double r1     = scalling * xx[0] / (1 - xx[0]);
  const double theta1 = 2 * M_PI * xx[1];
  const double phi1   = M_PI * xx[2];
  const std::vector<double> p1 =
      r1 * (std::vector<double>){
               cos(theta1) * sin(phi1), sin(theta1) * sin(phi1), cos(phi1)};

  // Spherical coordinates p2
  const double r2     = scalling * xx[3] / (1 - xx[3]);
  const double theta2 = 2 * M_PI * xx[4];
  const double phi2   = M_PI * xx[5];
  const std::vector<double> p2 =
      r2 * (std::vector<double>){
               cos(theta2) * sin(phi2), sin(theta2) * sin(phi2), cos(phi2)};

  const double E1 = proc->Energy(proc->m1, p1);
  const double E2 = proc->Energy(proc->m2, p2);

  if (r1 == 0 or r2 == 0 or phi1 == 0 or phi1 == M_PI or phi2 == 0 or
      phi2 == M_PI)
  {
    ff[0] = 0;
    return 0;
  }

  ff[0] = proc->MonteCarloInt(E1, E2, p1, p2) * _4_M_4 * r1 *
          pow(r1 + scalling, 2) * sin(phi1) * sin(phi2) * r2 *
          pow(r2 + scalling, 2) / (4 * scalling * scalling);

  if (proc->m1 > 0) ff[0] *= r1 / E1;
  if (proc->m2 > 0) ff[0] *= r2 / E2;

  return 0;
};
