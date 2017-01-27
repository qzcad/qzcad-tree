#ifndef RFUNCTIONS_H
#define RFUNCTIONS_H

namespace msh {
/**
 * @brief Коънюкция
 * @param x Первый аргумент
 * @param y Второй агрумент
 * @return x + y - sqrt(x*x + y*y)
 */
double con(const double &x, const double &y);
/**
 * @brief Дизъюнкция
 * @param x Первый аргумент
 * @param y Второй аргумент
 * @return x + y + sqrt(x*x + y*y)
 */
double dis(const double &x, const double &y);
/**
 * @brief Функция области, ограниченной окружностью с центром в начале координат
 * @param x Абсцисса
 * @param y Ордината
 * @param r Радиус
 * @return r^2 - x^2 - y^2
 */
double circle(const double &x, const double &y, const double &r);
/**
 * @brief Функция области, ограниченной эллипсом с центром в начале координат
 * @param x Абсцисса
 * @param y Ордината
 * @param a Большая полуось
 * @param b Малая полуось
 * @return 1 - x^2 / a^2 - y^2 / b^2
 */
double ellipse(const double &x, const double &y, const double &a, const double &b);
/**
 * @brief Полоса, перепендикулярная оси, симметричная относительно начала отсчета
 * @param x Координата оси
 * @param w Ширина полосы
 * @return (w/2)^2 - x^2
 */
double band(const double &x, const double &w);
/**
 * @brief Функция области, ограниченной прямоугольником с центром в начале координат
 * @param x Абсцисса
 * @param y Ордината
 * @param w Ширина
 * @param h Высота
 * @return Значение con(band(x, w), band(y, h))
 * @see con, band
 */
double rectangle(const double &x, const double &y, const double &w, const double &h);
/**
 * @brief Функция области, ограниченной прямоугольником с центром в начале координат
 * @param x Абсцисса
 * @param y Ордината
 * @param w Ширина
 * @param h Высота
 * @param r Радиус скругления углов
 * @return Значение функции области
 */
double rectangle(const double &x, const double &y, const double &w, const double &h, const double &r);
/**
 * @brief Область находящаяся "выше" прямой, определенной точками (x1, y1), (x2, y2)
 * @param x Абсцисса
 * @param y Ордината
 * @param x1 Абсцисса первой точки прямой
 * @param y1 Ордината первой точки прямой
 * @param x2 Абсцисса второй точки прямой
 * @param y2 Ордината второй точки прямой
 * @return (y - y1) * (x2 - x1) - (x - x1) * (y2 - y1)
 */
double line(const double &x, const double &y, const double &x1, const double &y1, const double &x2, const double &y2);
/**
 * @brief Область, ограниченная правильным многоуголником, вписаным в окружность
 * @param x Абсцисса
 * @param y Ордината
 * @param r Радиус окружности
 * @param n Количество вершин в многоугольнике
 * @return Значение соответствующей функции
 */
double regular(const double &x, const double &y, const double &r, const int &n);
/**
 * @brief Область, ограниченная эллипсоидом
 * @param x Абсцисса
 * @param y Ордината
 * @param z Аппликата
 * @param a Первая полуось
 * @param b Вторая полуось
 * @param c Третья полуось
 * @return Значение соответствующей функции
 */
double ellipsoid(const double &x, const double &y, const double &z, const double &a, const double &b, const double &c);
/**
 * @brief Область, ограниченная сферой
 * @param x Абсцисса
 * @param y Ордината
 * @param z Аппликата
 * @param r Радиус
 * @return Значение соответствующей функции
 */
double sphere(const double &x, const double &y, const double &z, const double &r);
/**
 * @brief Область, огранинченная плоскостью (ориентированное уравнение плоскости), заданной тремя точками
 * @param x Абсцисса
 * @param y Ордината
 * @param z Аппликата
 * @param x1 Абсцисса первой точки
 * @param y1 Ордината первой точки
 * @param z1 Аппликата первой точки
 * @param x2 Абсцисса второй точки
 * @param y2 Ордината второй точки
 * @param z2 Аппликата второй точки
 * @param x3 Абсцисса третей точки
 * @param y3 Ордината третей точки
 * @param z3 Аппликата третей точки
 * @return Значение соответствующей функции
 */
double plane(const double &x, const double &y, const double &z,
             const double &x1, const double &y1, const double &z1,
             const double &x2, const double &y2, const double &z2,
             const double &x3, const double &y3, const double &z3);
}

#endif // RFUNCTIONS_H
