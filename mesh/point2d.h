/**
  * @author Сергей Чопоров
  * @date 15/10/2012
  * @version 1.0.4
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
     * @brief Конструктор
     * @param x Ордината
     * @param y Абсцисаа
     */
    Point2D(double x = 0.0, double y = 0.0);
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
     * @brief Абсцисса
     * @return Абсциссу точки
     */
    virtual double y() const;
    /**
     * @brief Длина вектора
     * @return Длину вектора, определенного точкой
     */
    virtual double length() const;
    /**
     * @brief Установить координаты точки
     * @param x Ордината
     * @param y Абсцисса
     */
    void set(double x, double y);
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
private:
    double y_; //!< Абсцисса точки
};
}

#endif // POINT2D_H
