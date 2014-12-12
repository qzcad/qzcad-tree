/**
  * @author Сергей Чопоров
  * @date 15/10/2012
  * @version 1.0.3
  * @copyright Copyright 2012 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef POINT1D_H
#define POINT1D_H

#include "point.h"

namespace msh
{
/**
 * @brief Класс Point1D - точка в 1D
 * @see Point
 */
class Point1D: public Point
{
public:
    /**
     * @brief Конструктор по умолчанию инициализирует абсциссу нулем
     */
    Point1D();
    /**
     * @brief Конструктор
     * @param x Начальное значение ординаты
     */
    Point1D(const double &x);
    /**
     * @brief Конструктор копирования
     * @param point Объект-точка для копирования
     */
    Point1D(const Point1D &point);
    /**
     * @brief Вектор, заданный двумя точками (вершина вектора - разность координат второй и первой точек)
     * @param firstPoint Первая точка
     * @param secondPoint Вторая точка
     */
    Point1D(const Point1D &firstPoint, const Point1D &secondPoint);
    /**
     * @brief Размерность пространства
     * @return 1
     */
    virtual int dimension() const;
    /**
     * @brief Абсцисса
     * @return Абсциссу точки или ее проекции
     */
    virtual double x() const;
    /**
     * @brief Ордината (заглушка)
     * @return 0.0
     */
    virtual double y() const;
    /**
     * @brief Аппликата (заглушка)
     * @return 0.0
     */
    virtual double z() const;
    /**
     * @brief Длина
     * @return Длину вектора (отрезка)
     */
    virtual double length() const;
    /**
     * @brief Установить координаты точки
     * @param x Новое значение координаты точки
     */
    void set(const double &x);
    /**
     * @brief setX Установить значение абсциссы точки
     * @param x Новое значение абсциссы точки
     */
    void setX(const double &x);
    /// @name Операторы
    /// @{
    /// Присваивание
    Point1D &operator =(const Point1D& point);
    /// Проверка на равенство
    friend bool operator ==(const Point1D &leftPoint, const Point1D &rightPoint);
    /// Унарный минус
    friend Point1D operator -(const Point1D &point);
    /// Бинарный минус векторов, определенных точками leftPoint и rightPoint
    friend Point1D operator -(const Point1D &leftPoint, const Point1D &rightPoint);
    /// Бинарный плюс векторов, определенных точками LeftPoint и RightPoint
    friend Point1D operator +(const Point1D &leftPoint, const Point1D &rightPoint);
    /// Скалярное произведение векторов, определенных точками this и point
    double operator *(const Point1D &point) const;
    /// Произведение числа d на вектор, заданный точкой P
    friend Point1D operator *(double d, const Point1D &point);
    /// Отношение вектора, занного точкой point, к числу dec
    friend Point1D operator /(const Point1D &point, double dec);
    /// @}
    /**
     * @brief Координаты нормализованного вектора
     * @return Координаты нормализованного вектора (единичного вектора того же направления)
     */
    Point1D normalized() const;
    /**
     * @brief Расстояние до точки
     * @param point Точка, до которой необходимо определить расстояние
     * @return Расстояние до точки point
     */
    double distanceTo(const Point1D &point) const;
    /**
     * @brief Проверка на приблизительное равенство
     * @param point Точка, с которой производится сравнение
     * @param epsilon Точность
     * @return true, если точки приблизительно равны; false - иначе
     */
    bool isEqualTo(const Point1D &point, double epsilon = 0.000001) const;
private:
    double x_; //!< Абсцисса точки
};
}

#endif // POINT1D_H
