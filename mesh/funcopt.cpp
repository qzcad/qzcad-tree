#include "funcopt.h"

#include <iostream>
#include <math.h>

namespace msh {

std::vector<double> conjugateGradient(CoordinateFunction Functional, const std::vector<double> &x0, const double &h, const double &epsilon, int maxIter, bool messages)
{
    std::vector<double> xk = x0;
    std::vector<double> dxk = nabla(Functional, xk, h);
    std::vector<double> sk = dxk;
    double alpha = goldenRatio(Functional, -1.0, 1.0, xk, sk, epsilon, maxIter);
    for(std::vector<double>::iterator ixk = xk.begin(), isk = sk.begin(); ixk != xk.end(); ++ixk, ++isk)
    {
        (*ixk) = (*ixk) + alpha * (*isk);
    }
//    double f0 = Functional(xk);

    for (int k = 0; k < maxIter; k++)
    {
        std::vector<double> dxk1 = nabla(Functional, xk, h);
        double beta = 0.0, sum_old = 0.0, sum_new = 0.0;
        for(std::vector<double>::iterator idxk = dxk.begin(), idxk1 = dxk1.begin(); idxk != dxk.end(); ++idxk, ++idxk1)
        {
            sum_old += (*idxk) * (*idxk);
            sum_new += (*idxk1) * (*idxk1);
        }
        if (sum_old != 0.0) beta = sum_new / sum_old;
        for(std::vector<double>::iterator isk = sk.begin(), idxk1 = dxk1.begin(); isk != sk.end(); ++isk, ++idxk1)
        {
            (*isk) = (*idxk1) + beta * (*isk);
        }
        alpha = goldenRatio(Functional, -1.0, 1.0, xk, sk, epsilon, maxIter);
        double resid = 0.0;
        for(std::vector<double>::iterator ixk = xk.begin(), isk = sk.begin(); ixk != xk.end(); ++ixk, ++isk)
        {
            double ask = alpha * (*isk);
            (*ixk) = (*ixk) + ask;
            resid += ask * ask;
        }
        resid = sqrt(resid);
//        double f1 = Functional(xk);
//        double resid = fabs(f1 - f0);

        if (messages && k % 100 == 0) std::cout << "residual: " << resid << std::endl;
        if (resid < epsilon) break;

        dxk = dxk1;
//        f0 = f1;
    }
    return xk;
}

std::vector<double> nabla(CoordinateFunction Functional, const std::vector<double> &x, double h)
{
    std::vector<double> gradient(x.size(), 0.0);
    std::vector<double> var = x;
    double x_plus_h; // x + h
    double x_minus_h; // x - h
    double x_plus_h_h; // x + 2 * h
    double x_minus_h_h; // x - 2 * h
    double f_x_plus_h; // f(x + h)
    double f_x_minus_h; // f(x - h)
    double f_x_plus_h_h; // f(x + 2 * h)
    double f_x_minus_h_h; // f(x - 2 * h)
    double x_i; // temp value for x[i]

    for (std::vector<double>::size_type i = 0; i < var.size(); i++)
    {
        x_i = var[i];
        x_plus_h = x_i + h;
        x_minus_h = x_i - h;
        x_plus_h_h = x_plus_h + h;
        x_minus_h_h = x_minus_h - h;
        var[i] = x_plus_h;
        f_x_plus_h = Functional(var);
        var[i] = x_minus_h;
        f_x_minus_h = Functional(var);
        var[i] = x_plus_h_h;
        f_x_plus_h_h = Functional(var);
        var[i] = x_minus_h_h;
        f_x_minus_h_h = Functional(var);
        gradient[i] = (f_x_minus_h_h - 8.0 * f_x_minus_h + 8.0 * f_x_plus_h - f_x_plus_h_h) / (12.0 * h);
        var[i] = x_i;
    }
    return gradient;
}

double goldenRatio(CoordinateFunction Functional, const double &a, const double &b, const std::vector<double> &x0, const std::vector<double> &s, const double &epsilon, int maxIter)
{
    const double phi = (sqrt(5.0) + 1.0) / 2.0; // golden ratio
    double left = a;
    double right = b;
    double c = right - (right - left) / phi;
    double d = left + (right - left) / phi;
    int i = 0;

    while (i < maxIter && fabs(d - c) > epsilon)
    {
        if (lambda(Functional, x0, s, c) < lambda(Functional, x0, s, d))
        {
            right = d;
        }
        else
        {
            left = c;
        }
        c = right - (right - left) / phi;
        d = left + (right - left) / phi;
    }
    return (left + right) / 2.0;
}

double lambda(CoordinateFunction Functional, const std::vector<double> &x, const std::vector<double> &s, const double &lambda_val)
{
    std::vector<double> x_lambda(x.size());
    std::vector<double>::const_iterator ix, is;
    std::vector<double>::iterator il;
    for (ix = x.begin(), is = s.begin(), il = x_lambda.begin();
         ix != x.end() && is != s.end() && il != x_lambda.end();
         ++ix, ++is, ++il)
    {
        (*il) = (*ix) + lambda_val * (*is);
    }
    return Functional(x_lambda);
}

double norm2(const std::vector<double> &x)
{
    double norm = 0.0;
    for (auto v: x)
    {
        norm += v * v;
    }
    return sqrt(norm);
}

std::vector<double> descentGradient(CoordinateFunction Functional, const std::vector<double> &x0, const double &h, const double &epsilon, int maxIter, bool messages)
{
    std::vector<double> xk = x0;
    double fxk = Functional(xk);

    if (messages == true)
        std::cout << "Gradient Descent Method analysis" << std::endl;

    for(int k = 0; k < maxIter; k++)
    {
        std::vector<double> sk = nabla(Functional, xk, h);
        double lambda = goldenRatio(Functional, -1.0, 1.0, xk, sk, epsilon, maxIter);

        std::vector<double>::iterator ixk;
        std::vector<double>::iterator isk;
        for (ixk = xk.begin(), isk = sk.begin(); ixk != xk.end() && isk != sk.end(); ++ixk, ++isk)
            (*ixk) += lambda * (*isk);

        double fxk1 = Functional(xk);
        double residual = fabs(fxk1 - fxk);

        if(residual < epsilon)
            break;

        fxk = fxk1;

        if (messages && (k%100 == 0)) std::cout << "residual: " << residual << std::endl;
    }
    return xk;
}

}
