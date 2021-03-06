/**
  * @author Сергей Чопоров
  * @date 15/10/2012
  * @version 1.0.4
  * @copyright Copyright 2012 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef POINT2D_H
#define POINT2D_H

#include "point1d.h"

namespace msh
{
/**
 * @brief Класс Point2D - точка на плоскости (2D)
 * @see Point1D
 */
class Point2D: public Point1D
{
public:
    /**
     * @brief Конструктор по умолчанию иницализирует точку началом координат (0; 0)
     */
    Point2D();
    /**
     * @brief Конструктор
     * @param x Абсцисса
     * @param y Ордината
     */
    Point2D(const double &x, const double &y);
    /**
     * @brief Конструктор копирования
     * @param point Объект-точка для копирования
     */
    Point2D(const Point2D &point);
    /**
     * @brief Конструктор, который формирует точку как вершину вектора, задоного двумя точками
     * @param firstPoint Первая точка
     * @param secondPoint Вторая точка
     */
    Point2D(const Point2D &firstPoint, const Point2D &secondPoint);
    /**
     * @brief Размерность пространства
     * @return
     */
    virtual int dimension() const;
    /**
     * @brief Ордината точки
     * @return Значение ординаты точки
     */
    virtual double y() const;
    /**
     * @brief Длина вектора
     * @return Длину вектора, определенного точкой
     */
    virtual double length() const;
    /**
     * @brief Установить координаты точки
     * @param x Абсцисса
     * @param y Ордината
     */
    void set(const double &x, const double &y);
    /**
     * @brief Напечатать на стандартную консоль координаты
     */
    void print() const;
    /**
     * @brief Напечатать на стандартную консоль координаты с переходом на новую строку
     */
    void println() const;
    /**
     * @brief Установить значение ординаты
     * @param y Новое значение ординаты
     */
    void setY(const double &y);
    /// @name Операторы
    /// @{
    /// Присваивание
    Point2D &operator =(const Point2D &point);
    /// Проверка на равенство leftPoint и rightPont
    friend bool operator ==(const Point2D &leftPoint, const Point2D &rightPoint);
    /// Унарный минус
    friend const Point2D operator -(const Point2D &point);
    /// Бинарный минус векторов, определенных точками leftPoint и rightPoint
    friend const Point2D operator -(const Point2D &leftPoint, const Point2D &rightPoint);
    /// Бинарный плюс векторов, определенных точками leftPoint и rightPoint
    friend const Point2D operator +(const Point2D &leftPoint, const Point2D &rightPoint);
    /// Скалярное произведение векторов, определенных точками this и point
    double operator *(const Point2D &point) const;
    /// Произведение числа dec на вектор, заданный точкой point
    friend const Point2D operator *(double dec, const Point2D &point);
    /// Отношение вектора, занного точкой point, к числу dec
    friend const Point2D operator /(const Point2D &point, double dec);
    /**
     * @brief Оператор сравнения "меньше" с использованием лексикографического порядка
     * @param leftPoint Левая точка
     * @param rightPoint Правая точка
     * @return leftPoint < rightPoint
     */
    friend bool operator <(const Point2D &leftPoint, const Point2D &rightPoint);
    /// @}
    /**
     * @brief Координаты нормализованного вектора
     * @return Координаты нормализованного вектора (единичного вектора того же направления)
     */
    Point2D normalized() const;
    /**
     * @brief Расстояние до точки
     * @param point Точка, до которой необходимо определить расстояние
     * @return Расстояние до точки point
     */
    double distanceTo(const Point2D &point) const;
    /**
     * @brief Метод для вычисления расстояние от точки до отрезка
     * @param segment0 Начало отрезка
     * @param segment1 Конец отрезка
     * @return Расстояние до ближайшей точки отрезка
     */
    double distanceTo(const Point2D &segment0, const Point2D &segment1) const;
    /**
     * @brief Проверка на приблизительное равенство
     * @param point Точка, с которой производится сравнение
     * @param epsilon Точность
     * @return true, если точки приблизительно равны; false - иначе
     */
    bool isEqualTo(const Point2D &point, double epsilon = 0.00001) const;
    /**
     * @brief Векторное произведение векторов, определенных текущей точкой и заданной
     * @param point Координаты ветора, с которым необходимо найти векторное произведение
     * @return Значение векторного произведения в третьем направлении (на плоскости два других обнуляются)
     */
    double product(const Point2D &point) const;
    /**
     * @brief Относительное положение точки и отрезка
     */
    enum PointToSegment{
        LEFT,       //!< Слева
        RIGHT,      //!< Справа
        BEYOND,     //!< Впереди
        BEHIND,     //!< Позади
        BETWEEN,    //!< Между началом и концом
        ORIGIN,     //!< Начало
        DESTINATION //!< Конец
    };
    PointToSegment classify(const Point2D &p0, const Point2D &p1);
    /**
     * @brief Метод возвращает вектор, перепендикулярный текущему (x, y) _|_ (-y, x)
     * @return Экзмепляр перпендикулярного вектора
     */
    Point2D perpendicular() const;
    /**
     * @brief Масштабировать координаты точки на коэффициент
     * @param d Коэффициент масштабирования
     */
    virtual void scale(const double &d);
    /**
     * @brief Проверка пересекаются ли прямые P0P1 и Q0Q1
     * @param P0 Первая точка прямой P
     * @param P1 Вторая точка прямой P
     * @param Q0 Первая точка прямой Q
     * @param Q1 Вторая точка прямой Q
     * @param p Коэффицент, соответсвующий точке пресечения на прямой P ( X = P0 + p * (P1 - P0) )
     * @param q Коэффицент, соответсвующий точке пресечения на прямой q ( X = Q0 + q * (Q1 - Q0) )
     * @return
     */
    friend bool isCrossed(const Point2D &P0, const Point2D &P1, const Point2D &Q0, const Point2D &Q1, double &p, double &q);
    /**
     * @brief Вычислить угол B-this-C
     * @param B Первая точка угла
     * @param C Вторая точка угла
     * @return Угол B-this-C в радианах
     */
    double angle(const Point2D &B, const Point2D &C) const;
private:
    double y_; //!< Ордината точки
};
}

#endif // POINT2D_H
