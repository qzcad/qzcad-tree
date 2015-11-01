#ifndef FUNCOPT_H
#define FUNCOPT_H

#include <vector>
#include <functional>

namespace msh {

typedef std::function<double(const std::vector<double> &)> CoordinateFunction;

std::vector<double> conjugateGradient(CoordinateFunction Functional,
                                      const std::vector<double> &x0,
                                      const double &h,
                                      const double &epsilon, int maxIter = 300);

std::vector<double> descentGradient(CoordinateFunction Functional,
                                    const std::vector<double> &x0,
                                    const double &h,
                                    const double &epsilon, int maxIter = 1000);

std::vector<double> nabla(CoordinateFunction Functional, const std::vector<double> &x, double h);

double goldenRatio(CoordinateFunction Functional, const double &a, const double &b, const std::vector<double> &x0, const std::vector<double> &s, const double &epsilon, int maxIter = 300);

double lambda(CoordinateFunction Functional, const std::vector<double> &x, const std::vector<double> &s, const double &lambda_val);

double norm2(const std::vector<double> &x);
}

#endif // FUNCOPT_H
