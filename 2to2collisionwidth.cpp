#include "2to2collisionwidth.hpp"
#include <fstream>
#include <iostream>
using namespace TwoToTwoCollisionWidth;

void my_gsl_error_handler(const char *reason,
                          const char *file,
                          int line,
                          int gsl_errno)
{
  // std::cerr << "GSL Error: " << reason << " in " << file << ":" << line
  //           << " (Error Code: " << gsl_errno << ")" << std::endl;

  return;
}

double Process::integrate_phi_s(const double &theta)
{
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(100000);

  gsl_set_error_handler(&my_gsl_error_handler);

  auto ptr = [=](const double &phi) -> double
  { return integrand_theta_phi_s(theta, phi); };

  gsl_function_pp<decltype(ptr)> Fp(ptr);
  const gsl_function *F = static_cast<gsl_function *>(&Fp);

  double result, error;
  if (gsl_integration_qag(F,
                          0,
                          M_PI,
                          1e-4,
                          1e-4,
                          100000,
                          GSL_INTEG_GAUSS41,
                          w,
                          &result,
                          &error) != 0)

  {
    /*
    std::cout << "-- Integrate phi --\n";

    std::ofstream MyFile("filename.tsv");

    for (double p = 0; p < M_PI; p += 1e-5)
    {
      // std::cout << p << "\t" << ptr(p) << "\n";
      MyFile << std::setprecision(10) << p << "\t" << ptr(p) << "\n";
    }
    // std::cout << M_PI << "\t" << ptr(M_PI) << "\n";
    MyFile << M_PI << "\t" << ptr(M_PI) << "\n";

    MyFile.close();

    std::cout << "\n\nresults\t" << result << "\n";
    std::cout << "error\t" << error << "\n";
    std::cout << "rel\t" << error / result << "\n";
    std::cout << "theta\t" << theta << "\n";
    std::cout << "theta - pi\t" << std::setprecision(20) << theta - M_PI
              << "\n";
    std::cout << "E1\t" << E1_ << "\n";
    std::cout << "E2\t" << ET_ - E1_ << "\n";
    std::cout << "p1\t" << p1_[0] << "\t" << p1_[1] << "\t" << p1_[2] << "\t"
              << "\n";
    std::cout << "p2\t" << (p1p2_ - p1_)[0] << "\t" << (p1p2_ - p1_)[1] << "\t"
              << (p1p2_ - p1_)[2] << "\t"
              << "\n\n\n";

    std::cout << "\nE1 + E2\t" << ET_ << "\n";
    std::cout << "Minimum Energy\t" << (ET_ - Energy(m3 + m4, p1p2_)) << "\n";
    // exit(0);*/
  }
  gsl_integration_workspace_free(w);
  return result;
}

double Process::integrate_theta_t(const double &r)
{
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(100000);

  gsl_set_error_handler(&my_gsl_error_handler);

  auto ptr = [=](const double &theta) -> double
  { return integrand_r_theta_t(r, theta); };

  gsl_function_pp<decltype(ptr)> Fp(ptr);
  const gsl_function *F = static_cast<gsl_function *>(&Fp);

  double result, error;
  if (gsl_integration_qag(F,
                          0,
                          2 * M_PI,
                          1e-4,
                          1e-4,
                          100000,
                          GSL_INTEG_GAUSS41,
                          w,
                          &result,
                          &error) != 0)

  {
    /*
    std::cout << "-- Integrate theta --\n";

    std::ofstream MyFile("filename.tsv");

    for (double p = 0; p < M_PI; p += 1e-5)
    {
      // std::cout << p << "\t" << ptr(p) << "\n";
      MyFile << std::setprecision(10) << p << "\t" << ptr(p) << "\n";
    }
    // std::cout << M_PI << "\t" << ptr(M_PI) << "\n";
    MyFile << M_PI << "\t" << ptr(M_PI) << "\n";

    MyFile.close();

    std::cout << "\n\nresults\t" << result << "\n";
    std::cout << "error\t" << error << "\n";
    std::cout << "rel\t" << error / result << "\n";
    std::cout << "theta\t" << theta << "\n";
    std::cout << "theta - pi\t" << std::setprecision(20) << theta - M_PI
              << "\n";
    std::cout << "E1\t" << E1_ << "\n";
    std::cout << "E2\t" << ET_ - E1_ << "\n";
    std::cout << "p1\t" << p1_[0] << "\t" << p1_[1] << "\t" << p1_[2] << "\t"
              << "\n";
    std::cout << "p2\t" << (p1p2_ - p1_)[0] << "\t" << (p1p2_ - p1_)[1] << "\t"
              << (p1p2_ - p1_)[2] << "\t"
              << "\n\n\n";

    std::cout << "\nE1 + E2\t" << ET_ << "\n";
    std::cout << "Minimum Energy\t" << (ET_ - Energy(m3 + m4, p1p2_)) << "\n";
    // exit(0);*/
  }
  gsl_integration_workspace_free(w);
  return result;
}

double Process::Integrate_s(const double &E1,
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
  { return static_cast<Process *>(params)->integrate_phi_s(theta); };

  F.function = lambda;
  F.params   = this;

  double result, error;

  // gsl_integration_qags(&F, 0, 2 * M_PI, 1e-4, 1e-4, 10000, w, &result,
  // &error);

  gsl_integration_qag(&F,
                      0,
                      2 * M_PI,
                      1e-2,
                      1e-2,
                      10000,
                      GSL_INTEG_GAUSS21,
                      w,
                      &result,
                      &error);
  gsl_integration_workspace_free(w);
  return result;
}

double Process::Integrate_t(const double &E1,
                            const double &ED,
                            const std::vector<double> &p1,
                            const std::vector<double> &p3,
                            const std::vector<double> &p1p3)
{
  E1_   = E1;
  ED_   = ED;
  p1_   = p1;
  p3_   = p3;
  p1p3_ = p1p3;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);
  gsl_set_error_handler(&my_gsl_error_handler);
  gsl_function F;

  // Capture 'this' using another lambda
  auto lambda = [](double r, void *params) -> double
  { return static_cast<Process *>(params)->integrate_theta_t(r); };

  F.function = lambda;
  F.params   = this;

  double result, error;

  // gsl_integration_qags(&F, 0, 2 * M_PI, 1e-4, 1e-4, 10000, w, &result,
  // &error);

  gsl_integration_qagiu(&F, 0, 1e-2, 1e-2, 10000, w, &result, &error);
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

void Process::Calculate_a_b(const double &omega,
                            const double &k,
                            double &REa,
                            double &IMa,
                            double &REb,
                            double &IMb)
{
  return;
}

double Process::ReT1(const double &g,
                     const double &C,
                     const double &omega,
                     const double &k,
                     const double &mB,
                     const double &mF)
{
  // All possible poles
  std::vector<double> poles = {
      (-omega - k) / 2., (-omega + k) / 2., (omega - k) / 2., (omega + k) / 2.};
  // Sort poles
  std::sort(poles.begin(), poles.end());
  // Get only strictly positive
  std::vector<double> positive_poles;
  for (auto const &pole : poles)
    if (pole > 0) positive_poles.push_back(pole);

  // Add 0 at front and some slack at the end
  positive_poles.insert(positive_poles.begin(), 0.);
  positive_poles.push_back(positive_poles.back() * 2.);
  positive_poles.erase(unique(positive_poles.begin(), positive_poles.end()),
                       positive_poles.end());

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
          F, positive_poles.back(), 1e-6, 1e-6, 10000, w, &result, &error) != 0)
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
  r = result;

  gsl_integration_workspace *w1 = gsl_integration_workspace_alloc(10000);

  // convert std::vector to array
  size_t np = positive_poles.size();
  double points[np];
  std::copy(positive_poles.begin(), positive_poles.end(), points);

  double result1, error1;
  int gsl_status = gsl_integration_qagp(
      F, points, np, 1e-6, 1e-6, 10000, w1, &result1, &error1);
  if (gsl_status != 0)
  {

    if (error1 / result1 > 1e-3)
    {
      /*
      std::cout << "error\t\t\t\t" << error1 / result1 << "\n";
      std::cout << "GSL error: " << gsl_strerror(gsl_status) << "\n";
      std::cout << "Real part of T1 (2))\n";
      std::cout << "Poles\t";
      for (auto const &pole : points)
        std::cout << pole << "\t";
      std::cout << "\n";
      std::ofstream MyFile("ret1.tsv");
      for (auto p = positive_poles.front(); p < positive_poles.back();
           p += (positive_poles.back() - positive_poles.front()) / 10000.)
      {
        MyFile << std::setprecision(10) << p << "\t" << ptr(p) << "\n";
      }
      MyFile.close();
      std::cout << "result1\t" << result1 << "\n";
      std::cout << "error1\t" << error1 << "\n";
      std::cout << "rel\t" << error1 / result1 << "\n";
      std::cout << "T\t" << T << "\n";
      std::cout << "omega\t" << omega << "\n";
      std::cout << "k\t" << k << "\n";
      std::cout << "omega - k\t" << omega - k << "\n";
      std::cout << "r = " << r << "\n";
      // exit(0);
      //  return 0;*/
    }
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
  if (omega < 0) return ImT1(g, C, -omega, k, mB, mF);

  double r = 0;

  if (omega > k)
  {
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);

    gsl_set_error_handler(&my_gsl_error_handler);

    auto ptr = [=](const double &p) -> double
    {
      return (Distribution(mB, {p, 0, 0}, -1) * (p - omega) +
              p * Distribution(mF, {p, 0, 0}, 1));
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
      std::cout << "Imaginary part of T1 (1))\n";
      std::cout << "error\t" << error << "\n";
      std::cout << "omega\t" << omega << "\n";
      std::cout << "k\t" << k << "\n";
      std::cout << "omega - k\t" << omega - k << "\n";
      exit(0);
    }
    gsl_integration_workspace_free(w);
    r = result;
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
    r = result; // This sign comes from the retarded propagator

    auto ptr2 = [=](const double &p) -> double
    {
      return -(Distribution(mB, {p, 0, 0}, -1) * (p + omega) +
               p * Distribution(mF, {p, 0, 0}, 1));
    };
    gsl_function_pp<decltype(ptr2)> Fp2(ptr2);
    gsl_function *F2 = static_cast<gsl_function *>(&Fp2);

    double result2, error2;
    if (gsl_integration_qagiu(
            F2, (k - omega) / 2., 1e-6, 1e-6, 10000, w2, &result2, &error2) !=
        0)
    {
      std::cout << "Imaginary part of T1 (3))\n";
      std::cout << "error\t" << error2 << "\n";
      std::cout << "omega\t" << omega << "\n";
      std::cout << "k\t" << k << "\n";
      std::cout << "omega - k\t" << omega - k << "\n";
      exit(0);
    }
    gsl_integration_workspace_free(w2);
    r += result2;
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
  if (abs(omega) < 1e-4) return (g, C, 1e-3, k, mB, mF);
  // All possible poles
  std::vector<double> poles = {
      (-omega - k) / 2., (-omega + k) / 2., (omega - k) / 2., (omega + k) / 2.};
  // Sort poles
  std::sort(poles.begin(), poles.end());
  // Get only strictly positive
  std::vector<double> positive_poles;
  for (auto const &pole : poles)
    if (pole > 0) positive_poles.push_back(pole);

  // Add 0 at front and some slack at the end
  positive_poles.insert(positive_poles.begin(), 0.);
  positive_poles.push_back(positive_poles.back() * 2.);
  positive_poles.erase(unique(positive_poles.begin(), positive_poles.end()),
                       positive_poles.end());

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
          F, positive_poles.back(), 1e-6, 1e-6, 10000, w, &result, &error) != 0)
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
  r                             = result;
  gsl_integration_workspace *w1 = gsl_integration_workspace_alloc(10000);

  // convert std::vector to array
  size_t np = positive_poles.size();
  double points[np];
  std::copy(positive_poles.begin(), positive_poles.end(), points);

  double result1, error1;
  if (gsl_integration_qagp(
          F, points, np, 1e-4, 1e-4, 10000, w1, &result1, &error1) != 0)
  {
    /*std::cout << "Real part of T2 (2))\n";
    std::cout << "np\t" << np << "\n";
    for (int i = 0; i < np; i++)
    {
      std::cout << "Poles >" << i << "\t" << points[i] << "\n";
    }
    std::cout << "\n";

    for (double p = positive_poles.front(); p < positive_poles.back();
         p += (positive_poles.back() - positive_poles.front()) / 100.)
    {
      std::cout << p << "\t" << ptr(p) << "\t" << L1(p, omega, k) << "\n";
    }

    std::cout << "result1\t" << result1 << "\n";
    std::cout << "error1\t" << error1 << "\n";
    std::cout << "omega\t" << omega << "\n";
    std::cout << "k\t" << k << "\n";
    std::cout << "omega - k\t" << omega - k << "\n";*/
    //  exit(0);
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
  if (omega < 0) return -ImT2(g, C, -omega, k, mB, mF);

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);

  gsl_set_error_handler(&my_gsl_error_handler);

  auto ptr = [=](const double &p) -> double
  {
    return (Distribution(mF, {p, 0, 0}, 1) - Distribution(mB, {p, 0, 0}, -1));
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

  return g * g * C / (8. * M_PI * k) * (omega * omega - k * k) * result;
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

double Process::integrand_r_theta_t(const double &r, const double &theta)
{

  const double denominator =
      m1 * m1 + m3 * m3 + 2 * (p1_ * p3_ - Energy(m1, p1_) * Energy(m3, p3_));

  const double sqrtpart =
      pow(pow(m2, 2) - pow(m4, 2), 2) +
      4 * (pow(m2, 2) + pow(m4, 2) + 2 * pow(r, 2)) *
          (-pow(m1, 2) - pow(m3, 2) - p1_ * p3_ + E1_ * Energy(m3, p3_)) +
      4 * pow(-pow(m1, 2) - pow(m3, 2) - p1_ * p3_ + E1_ * Energy(m3, p3_), 2);

  const long double pzz1 =
      (abs(ED_) < 1)
          ? -0.5 *
                (-(pow(ED_, 2) * sqrt(pow(p1p3_[0], 2) + pow(p1p3_[2], 2))) +
                 sqrt(pow(p1p3_[0], 2) + pow(p1p3_[2], 2)) *
                     (-pow(m2, 2) + pow(m4, 2) + pow(p1p3_[0], 2) +
                      pow(p1p3_[2], 2)) +
                 ED_ *
                     sqrt(pow(ED_, 4) + pow(m2, 4) + pow(m4, 4) +
                          2 * pow(m4, 2) * pow(p1p3_[0], 2) + pow(p1p3_[0], 4) +
                          2 * pow(m4, 2) * pow(p1p3_[2], 2) +
                          2 * pow(p1p3_[0], 2) * pow(p1p3_[2], 2) +
                          pow(p1p3_[2], 4) +
                          2 * pow(m2, 2) *
                              (-pow(m4, 2) + pow(p1p3_[0], 2) +
                               pow(p1p3_[2], 2)) +
                          4 * pow(p1p3_[0], 2) * pow(r, 2) +
                          4 * pow(p1p3_[2], 2) * pow(r, 2) -
                          2 * pow(ED_, 2) *
                              (pow(m2, 2) + pow(m4, 2) + pow(p1p3_[0], 2) +
                               pow(p1p3_[2], 2) + 2 * pow(r, 2)))) /
                (pow(ED_, 2) - pow(p1p3_[0], 2) - pow(p1p3_[2], 2))
          : -0.5 *
                (-(pow(ED_, 2) * sqrt(pow(p1p3_[0], 2) + pow(p1p3_[2], 2))) +
                 sqrt(pow(p1p3_[0], 2) + pow(p1p3_[2], 2)) *
                     (-pow(m2, 2) + pow(m4, 2) + pow(p1p3_[0], 2) +
                      pow(p1p3_[2], 2)) +
                 ED_ * sqrt(sqrtpart)) /
                (pow(ED_, 2) - pow(p1p3_[0], 2) - pow(p1p3_[2], 2));
  const long double pzz2 =
      abs(ED_) < 1
          ? (pow(ED_, 2) * sqrt(pow(p1p3_[0], 2) + pow(p1p3_[2], 2)) -
             sqrt(pow(p1p3_[0], 2) + pow(p1p3_[2], 2)) *
                 (-pow(m2, 2) + pow(m4, 2) + pow(p1p3_[0], 2) +
                  pow(p1p3_[2], 2)) +
             ED_ *
                 sqrt(pow(ED_, 4) + pow(m2, 4) + pow(m4, 4) +
                      2 * pow(m4, 2) * pow(p1p3_[0], 2) + pow(p1p3_[0], 4) +
                      2 * pow(m4, 2) * pow(p1p3_[2], 2) +
                      2 * pow(p1p3_[0], 2) * pow(p1p3_[2], 2) +
                      pow(p1p3_[2], 4) +
                      2 * pow(m2, 2) *
                          (-pow(m4, 2) + pow(p1p3_[0], 2) + pow(p1p3_[2], 2)) +
                      4 * pow(p1p3_[0], 2) * pow(r, 2) +
                      4 * pow(p1p3_[2], 2) * pow(r, 2) -
                      2 * pow(ED_, 2) *
                          (pow(m2, 2) + pow(m4, 2) + pow(p1p3_[0], 2) +
                           pow(p1p3_[2], 2) + 2 * pow(r, 2)))) /
                (2. * (pow(ED_, 2) - pow(p1p3_[0], 2) - pow(p1p3_[2], 2)))
          : (pow(ED_, 2) * sqrt(pow(p1p3_[0], 2) + pow(p1p3_[2], 2)) -
             sqrt(pow(p1p3_[0], 2) + pow(p1p3_[2], 2)) *
                 (-pow(m2, 2) + pow(m4, 2) + pow(p1p3_[0], 2) +
                  pow(p1p3_[2], 2)) +
             ED_ * sqrt(sqrtpart)) /
                (2. * (pow(ED_, 2) - pow(p1p3_[0], 2) - pow(p1p3_[2], 2)));

  // const double pzz = (p1p3_[2] < 0) ? pzz1 : pzz2;
  const double pzz = pzz2;

  const double phi =
      acos(-p1p3_[2] / sqrt(pow(p1p3_[0], 2) + pow(p1p3_[2], 2)));

  const std::vector<double> p2 = {r * cos(phi) * cos(theta) + pzz * sin(phi),
                                  r * sin(theta),
                                  pzz * cos(phi) - r * cos(theta) * sin(phi)

  };

  const std::vector<double> p4 = p1p3_ + p2;

  const double delta_pzz =
      abs(pzz / sqrt(pow(m2, 2) + pow(pzz, 2) + pow(r, 2)) +
          (sqrt(pow(p1p3_[0], 2) + pow(p1p3_[2], 2)) - pzz) /
              sqrt(pow(m4, 2) + pow(p1p3_[0], 2) + pow(p1p3_[2], 2) -
                   2 * sqrt(pow(p1p3_[0], 2) + pow(p1p3_[2], 2)) * pzz +
                   pow(pzz, 2) + pow(r, 2)));

  if (abs((Energy(m1, p1_) + Energy(m2, p2)) /
              (Energy(m3, p3_) + Energy(m4, p4)) -
          1) > 1e-4)
  {
    std::cout << "------------------\tEnergy mismatch "
                 "(t-channel)\t---------------------\n";
    std::cout << std::setprecision(10) << "r\t" << r << "\n";
    std::cout << "theta\t" << theta << "\n\n";
    std::cout << "den\t" << denominator << "\t"
              << pow(ED_, 2) - pow(p1p3_[0], 2) - pow(p1p3_[2], 2) << "\n";
    std::cout << "sqrt part\t" << sqrtpart << "\n";
    std::cout << "phi\t" << phi << "\n";
    std::cout << "pzz1\t" << pzz1 << "\n";
    std::cout << "pzz2\t" << pzz2 << "\n";
    std::cout << "pzz\t" << pzz << "\n";
    std::cout << "ED\t" << ED_ << "\n";
    std::cout << "p1p3\t" << p1p3_[0] << "\t" << p1p3_[1] << "\t" << p1p3_[2]
              << "\n";
    std::cout << "delta_r\t" << delta_pzz << "\n";
    std::cout << "p1\t" << p1_[0] << "\t" << p1_[1] << "\t" << p1_[2] << "\n";
    std::cout << "p2\t" << p2[0] << "\t" << p2[1] << "\t" << p2[2] << "\n";
    std::cout << "p3\t" << p3_[0] << "\t" << p3_[1] << "\t" << p3_[2] << "\n";
    std::cout << "p4\t" << p4[0] << "\t" << p4[1] << "\t" << p4[2] << "\n";

    std::cout << "\n--------\n";

    std::cout << "E1\t" << E1_ << "\n";
    std::cout << "E2\t" << Energy(m2, p2) << "\n";
    std::cout << "E3\t" << E1_ - ED_ << "\n";
    std::cout << "E4\t" << Energy(m4, p4) << "\n";

    std::cout << "\nE1 + E2\t" << Energy(m1, p1_) + Energy(m2, p2) << "\n";
    std::cout << "E3 + E4\t" << Energy(m3, p3_) + Energy(m4, p4) << "\n";

    std::cout << "--------\n";

    std::cout << "(E1 + E2) / (E3 + E4) - 1\t"
              << abs((Energy(m1, p1_) + Energy(m2, p2)) /
                         (Energy(m3, p3_) + Energy(m4, p4)) -
                     1)
              << "\n";
    std::cout << "E1 + E2 - (E3 + E4)\t"
              << (Energy(m1, p1_) + Energy(m2, p2)) -
                     (Energy(m3, p3_) + Energy(m4, p4))
              << "\n";

    std::cout << "--------\n";

    std::cout << "m1\t" << m1 << "\n";
    std::cout << "m2\t" << m2 << "\n";
    std::cout << "m3\t" << m3 << "\n";
    std::cout << "m4\t" << m4 << "\n";

    std::cout << "--------\n";

    std::cout << "r\t" << r << "\n";
    std::cout << "theta\t" << theta << "\n";

    throw;
  }

  /*double mtinf = gs / sqrt(6.) * T; // Top thermal mass

  const double omega = Energy(0, p1_) - Energy(0, p3);
  const double k     = Energy(0, p1_ - p3);
  const double t =
      -2 * Energy(0, p1_) * Energy(0, p3) + 2 * p1_ * p3; // omega^2-k^2
  const double a = HTLa(mtinf, omega, k);
  const double b = HTLb(mtinf, omega, k);

  std::cout << "Phi " << phi << " |\t" << AmplitudeSquared_t(p1_, p2_, p3) <<
     "\t"
            << pow(-pow(mtinf, 2) + t, 2) << "\t"
            << pow(pow(b, 2) + 2 * (1 + a) * b * omega + pow(1 + a, 2) * t, 2)
            << "\n";*/

  if (isinf(AmplitudeSquared_t(p1_, p2, p3_) * Distribution(m2, p2, s2) *
            (1 + Distribution(m4, p4, s4)) * r /
            (delta_pzz * 4 * Energy(m2, p2) * Energy(m4, p4))))
  {
    std::cout << "\n-----\tAmpltidude is divergent\t-----\n\n";
    std::cout << "p1\t" << p1_[0] << "\t" << p1_[1] << "\t" << p1_[2] << "\n";
    std::cout << "p2\t" << p2[0] << "\t" << p2[1] << "\t" << p2[2] << "\n";
    std::cout << "p3\t" << p3_[0] << "\t" << p3_[1] << "\t" << p3_[2] << "\n";
    std::cout << "p4\t" << p4[0] << "\t" << p4[1] << "\t" << p4[2] << "\n\n";

    std::cout << "r\t" << r << "\n";
    std::cout << "theta\t" << theta << "\n";
    std::cout << "delta_pzz\t" << delta_pzz << "\n";

    std::cout << "Amp\t" << AmplitudeSquared_t(p1_, p2, p3_) << "\n";
    std::cout << "Distribution(m2, p2, s2)\t" << Distribution(m2, p2, s2)
              << "\n";
    std::cout << " Distribution(m4, p4, s4)\t" << Distribution(m4, p4, s4)
              << "\n";
    std::cout << "Amp\t" << AmplitudeSquared_t(p1_, p2, p3_) << "\n";
    std::cout << "Energy(m2, p2)\t" << Energy(m2, p2) << "\n";
    std::cout << "Energy(m4, p4)\t" << Energy(m4, p4) << "\n";
    exit(0);
  }

  return AmplitudeSquared_t(p1_, p2, p3_) * Distribution(m2, p2, s2) *
         (1 + Distribution(m4, p4, s4)) * r /
         (delta_pzz * 4 * Energy(m2, p2) * Energy(m4, p4));
}

double Process::integrand_theta_phi_s(const double &theta, const double &phi)
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

  if (abs(ET_ / (Energy(m3, p3) + Energy(m4, p4)) - 1) > 1e-6)
  {
    std::cout << "------------------\tEnergy mismatch "
                 "(s-channel)\t---------------------\n";
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

    std::cout << "(E1 + E2) / (E3 + E4) - 1\t"
              << abs(ET_ / (Energy(m3, p3) + Energy(m4, p4)) - 1) << "\n";
    std::cout << "E1 + E2 - (E3 + E4)\t"
              << ET_ - (Energy(m3, p3) + Energy(m4, p4)) << "\n";

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

  /*double mtinf = gs / sqrt(6.) * T; // Top thermal mass

  const double omega = Energy(0, p1_) - Energy(0, p3);
  const double k     = Energy(0, p1_ - p3);
  const double t =
      -2 * Energy(0, p1_) * Energy(0, p3) + 2 * p1_ * p3; // omega^2-k^2
  const double a = HTLa(mtinf, omega, k);
  const double b = HTLb(mtinf, omega, k);

  std::cout << "Phi " << phi << " |\t" << AmplitudeSquared_t(p1_, p2_, p3) <<
     "\t"
            << pow(-pow(mtinf, 2) + t, 2) << "\t"
            << pow(pow(b, 2) + 2 * (1 + a) * b * omega + pow(1 + a, 2) * t, 2)
            << "\n";*/

  return AmplitudeSquared_t(p1_, p2_, p3) * (1 - Distribution(m3, p3, s3)) *
         (1 + Distribution(m4, p4, s4)) * r * r * sin(phi) /
         (delta_r * 4 * Energy(m3, p3) * Energy(m4, p4));
}

double Process::MonteCarloInt_s(const double &E1,
                                const double &E2,
                                const std::vector<double> p1,
                                const std::vector<double> p2)
{

  const std::vector<double> p1p2 = p1 + p2;

  // Theta function to see if its possible to have a
  // process with a given p1 and p2
  if (E1 + E2 - Energy(m3 + m4, p1p2) <= 0) return 0.;

  double integral = Integrate_s(E1, E1 + E2, p1, p2, p1p2);

  // Normalization
  integral *=
      _2_PI_FACTORS * Distribution(m1, p1, s1) * Distribution(m2, p2, s2);

  return integral / prefactor;
};

double Process::MonteCarloInt_t(const double &E1,
                                const double &E3,
                                const std::vector<double> p1,
                                const std::vector<double> p3)
{
  // See if its possible to have a
  // process with a given p1 and p2
  const double p1p3 = sqrt((p1 - p3) * (p1 - p3));
  const double ED   = E1 - E3;

  //  There are three maximum value points. We need at least two different
  //  signs.
  double cond1 = ED + p1p3;
  double cond2 = ED - p1p3;
  double cond3 =
      (m2 == m4)
          ? 1
          : ED +
                sqrt(pow(m2, 2) +
                     pow((pow(m2, 2) * p1p3) / (pow(m2, 2) - pow(m4, 2)) +
                             (m2 * m4 * p1p3) / (pow(m2, 2) - pow(m4, 2)),
                         2)) -
                m4 * sqrt(1 + pow(-((m2 * p1p3) / (pow(m2, 2) - pow(m4, 2))) +
                                      (m4 * p1p3) / (-pow(m2, 2) + pow(m4, 2)),
                                  2)); // If m2 = m4 this does not matter

  if (cond1 * cond2 * cond3 == 0)
  {
    std::cout << "Error with zero cond1, cond2 or cond3\n";
    std::cout << "cond1 = " << cond1 << "\n";
    std::cout << "cond2 = " << cond2 << "\n";
    std::cout << "cond3 = " << cond3 << "\n";
  }

  if (cond1 * cond2 > 0 and cond1 * cond3 > 0)
    return 0.; // Then cond2 * cond3 > 0 and all have the same sign

  // Calculate a and b for the t-channel
  const double omega = Energy(0, p1) - Energy(0, p3);
  const double k     = Energy(0, p1 - p3);

  Calculate_a_b(omega, k, REa_t, IMa_t, REb_t, IMb_t);

  double integral = Integrate_t(E1, ED, p1, p3, p1 - p3);

  // std::cout << "integral is\t" << integral << "\n";

  if (isinf(integral))
  {
    std::cout << "p1\t" << p1[0] << "\t" << p1[1] << "\t" << p1[2] << "\n";
    std::cout << "p3\t" << p3[0] << "\t" << p3[1] << "\t" << p3[2] << "\n";
    std::cout << "yo?\n";
    exit(0);
  }

  // Normalization
  integral *=
      _2_PI_FACTORS * Distribution(m1, p1, s1) * (1 - Distribution(m3, p3, s3));

  return integral / prefactor;
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

int Process::Integrand_s(const int *ndim,
                         const cubareal xx[],
                         const int *ncomp,
                         cubareal ff[],
                         void *userdata)
{
  auto proc = static_cast<Process *>(userdata);

  if (xx[0] == 1 or xx[1] == 1)
  {
    // TODO might be removable
    ff[0] = 0.;
    return 0;
  }

  const double scalling =
      proc->T; // Put the hypercube on the region of interest.

  // Spherical coordinates p1
  const double r1              = scalling * xx[0] / (1 - xx[0]);
  const std::vector<double> p1 = (std::vector<double>){0, 0, r1};

  // Spherical coordinates p2
  const double r2  = scalling * xx[1] / (1 - xx[1]);
  const double phi = M_PI * xx[2];
  const std::vector<double> p2 =
      r2 * (std::vector<double>){sin(phi), 0., cos(phi)};

  const double E1 = proc->Energy(proc->m1, p1);
  const double E2 = proc->Energy(proc->m2, p2);

  if (r1 == 0 or r2 == 0 or phi == 0 or phi == M_PI)
  {
    ff[0] = 0;
    return 0;
  }

  ff[0] = proc->MonteCarloInt_s(E1, E2, p1, p2) * _4_M_4 * r1 *
          pow(r1 + scalling, 2) * sin(phi) * r2 * pow(r2 + scalling, 2) /
          (4 * scalling * scalling);

  if (proc->m1 > 0) ff[0] *= r1 / E1;
  if (proc->m2 > 0) ff[0] *= r2 / E2;

  return 0;
};

int Process::Integrand_t(const int *ndim,
                         const cubareal xx[],
                         const int *ncomp,
                         cubareal ff[],
                         void *userdata)
{
  auto proc = static_cast<Process *>(userdata);

  if (xx[0] == 1 or xx[1] == 1)
  {
    // TODO might be removable
    ff[0] = 0.;
    return 0;
  }

  const double scalling =
      proc->T; // Put the hypercube on the region of interest.

  // Spherical coordinates p1
  const double r1              = scalling * xx[0] / (1 - xx[0]);
  const std::vector<double> p1 = (std::vector<double>){0, 0, r1};

  // Spherical coordinates p3
  const double r3  = scalling * xx[1] / (1 - xx[1]);
  const double phi = M_PI * xx[2];
  const std::vector<double> p3 =
      r3 * (std::vector<double>){sin(phi), 0., cos(phi)};

  const double E1 = proc->Energy(proc->m1, p1);
  const double E3 = proc->Energy(proc->m3, p3);

  if (r1 == 0 or r3 == 0 or phi == 0 or phi == M_PI)
  {
    ff[0] = 0;
    return 0;
  }

  ff[0] = proc->MonteCarloInt_t(E1, E3, p1, p3) * _4_M_4 * r1 *
          pow(r1 + scalling, 2) * sin(phi) * r3 * pow(r3 + scalling, 2) /
          (4 * scalling * scalling);

  if (proc->m1 > 0) ff[0] *= r1 / E1;
  if (proc->m3 > 0) ff[0] *= r3 / E3;

  return 0;
};
