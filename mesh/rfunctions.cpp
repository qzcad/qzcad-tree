#include "rfunctions.h"

#include <math.h>

namespace msh {

double con(const double &x, const double &y)
{
    return x + y - sqrt(x*x + y*y);
}

double dis(const double &x, const double &y)
{
    return x + y + sqrt(x*x + y*y);
}

double circle(const double &x, const double &y, const double &r)
{
    return r*r - x*x - y*y;
}

double band(const double &x, const double &w)
{
    double w_2 = 0.5 * w;
    return w_2*w_2 - x*x;
}

double rectangle(const double &x, const double &y, const double &w, const double &h)
{
    return con(band(x, w), band(y, h));
}

double line(const double &x, const double &y, const double &x1, const double &y1, const double &x2, const double &y2)
{
    return (y - y1) * (x2 - x1) - (x - x1) * (y2 - y1);
}

double regular(const double &x, const double &y, const double &r, const int &n)
{
    double a = 2.0 * M_PI / (double)n;
    double res = line(x, y, r, 0.0, r * cos(a), r * sin(a));
    double alpha = a;
    double x1 = r * cos(a);
    double y1 = r * sin(a);
    double x2 = 0.0;
    double y2 = 0.0;
    for (int i = 1; i < n; i++)
    {
        alpha += a;
        x2 = r * cos(alpha);
        y2 = r * sin(alpha);
        res = con(res, line(x, y, x1, y1, x2, y2));
        x1 = x2;
        y1 = y2;
    }
    return res;
}

}
