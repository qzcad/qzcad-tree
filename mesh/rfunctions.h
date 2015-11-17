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
 * @brief Область находящаяся "ниже" прямой, определенной точками (x1, y1), (x2, y2)
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
}

#endif // RFUNCTIONS_H