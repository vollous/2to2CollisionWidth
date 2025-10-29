#include "numerics.hpp"

double Polylog(const int &s, const double &x)
{
  if (x > 1) throw("Polylog does not converge.");
  double r  = pow(x, s);
  double rr = r;
  for (int k = 2; k <= 10000; k++)
  {
    r += pow(x, k) / pow(k, s);
    if (abs(r / rr - 1) < 1e-7) return rr;
    rr = r;
  }
  std::cout << "Convergence issues on the polylog\n";
  return r;
}