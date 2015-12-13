#include <iostream>
#include <list>
#undef __STRICT_ANSI__
#include <math.h>
#include <stdarg.h>

#include <omp.h>

#include "trianglemesh2d.h"

using namespace std;
/**
 * @brief Функция для вычисления квадрата числа
 * @param a Число, квадрат которого необходимо вычислить
 * @return
 */
double square(double a)
{
    return a * a;
}
/**
 * @brief Функция вычисления конъюнкции
 * @param num Количество аргументов
 * @return Конъюнкцию аргументов
 */
double con( int num, ... )
{
    va_list arguments;
    double sum = 0.0;

    va_start ( arguments, num );
    sum = va_arg ( arguments, double );
    for ( int i = 1; i < num; i++ )
    {
        double x = va_arg ( arguments, double );
        sum += (x - sqrt(sum*sum + x*x));
    }
    va_end ( arguments );

    return sum;
}
/**
 * @brief Функция для вычисления дизъюнкции
 * @param n_args Количество аргументов
 * @return Дизъюнкция аргументов
 */
double dis( int num, ... )
{
    va_list arguments;
    double sum = 0.0;

    va_start ( arguments, num );
    sum = va_arg ( arguments, double );
    for ( int i = 1; i < num; i++ )
    {
        double x = va_arg ( arguments, double );
        sum += (x + sqrt(sum*sum + x*x));
    }
    va_end ( arguments );

    return sum;
}

double rect(double x, double y, double w, double h)
{
    double a = w / 2.0;
    double b = h / 2.0;
    return con(2, a*a - x*x, b*b - y*y);
}

double circle(double x, double y, double r)
{
    return r*r - x*x - y*y;
}

double planka1(double x, double y)
{
    double w = 30.0;
    double h = 20.0;
    double r = 5.0;
    return con(2, dis(2, rect(x, y, w, h), circle(x - w / 2.0, y, h / 2.0)), -circle(x, y, r));
}

msh::TriangleMesh2D *meshPlanka1(int count)
{
    double w = 30.0;
    double h = 20.0;
    list<msh::Point2D> character;
    character.push_back(msh::Point2D(-w / 2.0, -h / 2.0));
    character.push_back(msh::Point2D(w / 2.0, -h / 2.0));
    character.push_back(msh::Point2D(w / 2.0, h / 2.0));
    character.push_back(msh::Point2D(-w / 2.0, h / 2.0));
    double time = omp_get_wtime();
    msh::TriangleMesh2D *triangles = new msh::TriangleMesh2D();
    triangles->functionalDomain(count, count, -w, -w, 2.0 * w, 2.0 * w, planka1, character);
    //    double xmin = -w / 2.0 - h / 2;
    //    double ymin = -w / 2.0 - h / 2;
    //    double hh = (w + h) / (count - 1.0);
    //    double sum = 0.0;
    //    double time = omp_get_wtime();
    //#pragma omp parallel for reduction(+:sum)
    //    for (int i = 0; i < count; i++)
    //    {
    //        double x = xmin + i * hh;
    //        for (int j = 0; j < count; j++)
    //        {
    //            double y = ymin + j * hh;
    //            double s = planka1(x, y);
    //            sum += s;
    //        }
    //    }
    cout << endl << "h = " << 2.0 * w / (count - 1.0) << endl;
    cout << "Time: " << omp_get_wtime() - time << " s" << endl;
    return triangles;
}

double prokladka(double x, double y)
{
    double r1 = 40.0;
    double r2 = 20.0;
    double r3 = 30.0;
    double r4 = 5.0;
    double r5 = 24.0;
    double r6 = 14.0;
    double w1 = 14.0;
    double h1 = 5.0;
    double w2 = 2.0 * r5;
    double h2 = 100.0;
    double h3 = 48.0;
    double c1x = r3 * cos(-M_PI / 4.0);
    double c1y = r3 * sin(-M_PI / 4.0);
    double c2x = r3 * cos(M_PI / 4.0);
    double c2y = r3 * sin(M_PI / 4.0);
    double c3x = r3 * cos(3.0 * M_PI / 4.0);
    double c3y = r3 * sin(3.0 * M_PI / 4.0);
    double c4x = r3 * cos(5.0 * M_PI / 4.0);
    double c4y = r3 * sin(5.0 * M_PI / 4.0);

    double cc = dis(3, circle(x, y, r1), circle(x, y + h2, r5), rect(x, y + h2 / 2.0, w2, h2));
    double ccc = con(10, cc,
                     -circle(x, y, r2),
                     -rect(x, y - r2, w1, 2.0 * h1),
                     -circle(x - c1x, y - c1y, r4),
                     -circle(x - c2x, y - c2y, r4),
                     -circle(x - c3x, y - c3y, r4),
                     -circle(x - c4x, y - c4y, r4),
                     -rect(x, y + h2 - h3 / 2.0, 2.0 * r6, h3),
                     -circle(x, y + h2, r6),
                     -circle(x, y + h2 - h3, r6));
    return ccc;
}

msh::TriangleMesh2D *meshPlrokladka(int count)
{
    double r1 = 40.0;
    double r2 = 20.0;
    double r5 = 24.0;
    double w1 = 14.0;
    double h1 = 5.0;
    double h2 = 100.0;
    list<msh::Point2D> character;
    character.push_back(msh::Point2D(-w1 / 2.0, sqrt(r2*r2 - w1*w1 / 4.0)));
    character.push_back(msh::Point2D( w1 / 2.0, sqrt(r2*r2 - w1*w1 / 4.0)));
    character.push_back(msh::Point2D(-w1 / 2.0, r2 + h1));
    character.push_back(msh::Point2D( w1 / 2.0, r2 + h1));
    double time = omp_get_wtime();
    msh::TriangleMesh2D *triangles = new msh::TriangleMesh2D();
    triangles->functionalDomain(count, 2 * count + 1, -r1, -h2 - r5, 2.0 * r1, r1 + h2 + r5, prokladka, character);
//    cout << endl << "hx = " << 2.0 * r1 / (count - 1.0) << " hy = " << (r1 + h2 + r5) / (2.0 * count) << endl;
    cout << "Time: " << omp_get_wtime() - time << " s" << endl;
    return triangles;
}


double perfor(double x, double y)
{
    double w = 80.0;
    double r = 5.0;
    double re = rect(x, y, w, w);
    double c = 0.0;
    double d[] = { -30.0, -15.0, 0.0, 15.0, 30.0 };
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            if (i == 0 && j == 0)
            {
                c = circle(x - d[i], y - d[j], r);
            } else
            {
                c = dis(2, c, circle(x - d[i], y - d[j], r));
            }
        }
    }
    return con(2, re, -c);
}

msh::TriangleMesh2D *meshPerfor(int count)
{
    double w = 80.0;
    list<msh::Point2D> character;
    character.push_back(msh::Point2D(-w / 2.0, -w / 2.0));
    character.push_back(msh::Point2D(-w / 2.0, w / 2.0));
    character.push_back(msh::Point2D(w / 2.0, w / 2.0));
    character.push_back(msh::Point2D(w / 2.0, -w / 2.0));
    double time = omp_get_wtime();
    msh::TriangleMesh2D *triangles = new msh::TriangleMesh2D();
    triangles->functionalDomain(count, count, -w / 2.0, -w / 2.0, w, w, perfor, character);
//    cout << endl << "hx = " << 2.0 * r1 / (count - 1.0) << " hy = " << (r1 + h2 + r5) / (2.0 * count) << endl;
    cout << "Time: " << omp_get_wtime() - time << " s" << endl;
    return triangles;
}

int main()
{
    omp_set_num_threads(4);

    for (int c = 21; c <= 101; c += 20)
    {
//        msh::TriangleMesh2D *triangles = meshPerfor(c);
        msh::TriangleMesh2D *triangles = meshPlanka1(c);
        delete triangles;
    }
    return 0;
}

