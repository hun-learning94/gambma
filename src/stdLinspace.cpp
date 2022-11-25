#include <vector>
#include "stdLinspace.h"

std::vector<double> std_linspace(double a, double b, std::size_t N)
{
  double h = (b - a) / static_cast<double>(N-1);
  std::vector<double> xs(N);
  std::vector<double>::iterator x;
  double val;
  for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h) {
    *x = val;
  }
  return xs;
}