#include "funcopt.h"

#include <iostream>
#include <math.h>

namespace msh {

std::vector<double> conjugateGradient(CoordinateFunction Functional, const std::vector<double> &x0, const double &h, const double &epsilon, int maxIter, bool messages)
{
    std::vector<double> xk = x0;
    std::vector<double>::size_type size = x0.size();

    for (int k = 0; k < maxIter; k++)
    {
        std::vector<double> xkj = xk;
        std::vector<double> xkj1(size);
        std::vector<double> dxkj(size);
        if (messages) std::cout << "niter: " << k << " :: ";

        for (int j = 0; j < maxIter; j++)
        {
            std::vector<double> Skj = nabla(Functional, xkj, h);
            double lambda = goldenRatio(Functional, -1.0, 1.0, xkj, Skj, epsilon, maxIter);

            for(std::vector<double>::size_type i = 0; i < size; i++)
                xkj1[i] = xkj[i] + lambda * Skj[i];

            std::vector<double> nablaxkj = nabla(Functional, xkj, h);

            double nkj = norm2(xkj);
            double nkj1 = norm2(xkj1);
            double omega = (nkj1 * nkj1) / (nkj * nkj);
            double ndxkj = 0.0;

            for(std::vector<double>::size_type i = 0; i < size; i++)
            {
                Skj[i] = -nablaxkj[i] + omega * Skj[i];
                dxkj[i] = xkj1[i] - xkj[i];
                ndxkj += dxkj[i]*dxkj[i];
            }

            if(sqrt(ndxkj) < epsilon)
            {
                return xkj1;
            }

            xkj = xkj1;

            if (messages && j%10 == 0) std::cout << '*';
        }
        xk = xkj;
        if (messages) std::cout << "residual: " << norm2(dxkj) << std::endl;
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
    double left = a;
    double right = b;
    const double phi = (sqrt(5.0) + 1.0) / 2.0; // golden ratio
    double x1 = b - (b - a) / phi;
    double x2 = a + (b - a) / phi;
    double y1 = lambda(Functional, x0, s, x1);
    double y2 = lambda(Functional, x0, s, x2);

    for(int count = 0; count < maxIter; count++)
    {
        if(y1 <= y2)
        {
            right = x2;
            x2 = x1;
            x1 = right - (right - left) / phi;
            y2 = y1;
            y1 = lambda(Functional, x0, s, x1);
        }
        else
        {
            left = x1;
            x1 = x2;
            x2 = left + (right - left) / phi;
            y1 = y2;
            y2 = lambda(Functional, x0, s, x2);
        }
        if((right - left) < epsilon) break;
    }
    if(y1 <= y2) return x1;
    return x2;
}

double lambda(CoordinateFunction Functional, const std::vector<double> &x, const std::vector<double> &s, const double &lambda_val)
{
    std::vector<double> x_lambda(x.size());
    for (std::vector<double>::size_type i = 0; i < x.size(); i++)
    {
        x_lambda[i] = x[i] + lambda_val * s[i];
    }
    return Functional(x_lambda);
}

double norm2(const std::vector<double> &x)
{
    double norm = 0.0;
    for (auto it = x.cbegin(); it != x.cend(); ++it)
    {
        norm += (*it) * (*it);
    }
    return norm;
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

//        for(std::vector<double>::size_type i = 0; i < xk.size(); i++)
//            xk[i] += lambda * Sk[i];
        std::vector<double>::iterator ixk;
        std::vector<double>::iterator isk;
        for (ixk = xk.begin(), isk = sk.begin(); ixk != xk.end() && isk != sk.end(); ++ixk, ++isk)
            (*ixk) += lambda * (*isk);

        double fxk1 = Functional(xk);
        double residual = fabs(fxk1 - fxk);

        if(residual < epsilon)
            break;

        fxk = fxk1;

        if (messages && (k%10 == 0)) std::cout << "residual: " << residual << std::endl;
    }
    return xk;
}

}
