#include "rfunctions.h"
#undef __STRICT_ANSI__
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
    return 1.0 - (x*x + y*y) / (r*r);
}

double band(const double &x, const double &w)
{
    return 1.0 - 4.0 * x*x / (w*w);
}

double rectangle(const double &x, const double &y, const double &w, const double &h)
{
    return con(band(x, w), band(y, h));
}

double line(const double &x, const double &y, const double &x1, const double &y1, const double &x2, const double &y2)
{
    double dx = x2 - x1;
    double dy = y2 - y1;
    return ((y - y1) * dx - (x - x1) * dy) / sqrt(dx*dx + dy*dy); // нормальное уравнение прямой
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

double ellipse(const double &x, const double &y, const double &a, const double &b)
{
    return 1.0 - x*x / (a*a) - y*y / (b*b);
}

double rectangle(const double &x, const double &y, const double &w, const double &h, const double &r)
{
    double rr = dis(rectangle(x, y, w - 2.0 * r, h), rectangle(x, y, w, h - 2.0 * r));
    double cc = dis(dis(dis(circle(x - (w/2.0 - r), y - (h/2.0 - r), r),
                            circle(x + (w/2.0 - r), y - (h/2.0 - r), r)),
                        circle(x + (w/2.0 - r), y + (h/2.0 - r), r)),
                    circle(x - (w/2.0 - r), y + (h/2.0 - r), r));
    return dis(rr, cc);
}

double plane(const double &x, const double &y, const double &z, const double &x1, const double &y1, const double &z1, const double &x2, const double &y2, const double &z2, const double &x3, const double &y3, const double &z3)
{
    double a = (y2 - y1) * (z3 - z1) - (z2 - z1) * (y3 - y1);
    double b = (z2 - z1) * (x3 - x1) - (x2 - x1) * (z3 - z1);
    double c = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);
    return (a * (x - x1) + b * (y - y1) + c * (z - z1)) / sqrt(a*a + b*b + c*c);
}

double ellipsoid(const double &x, const double &y, const double &z, const double &a, const double &b, const double &c)
{
    return 1.0 - x*x / (a*a) - y*y / (b*b) - z*z / (c*c);
}

double sphere(const double &x, const double &y, const double &z, const double &r)
{
    return 1.0 - (x*x + y*y + z*z) / (r*r);
}

}
