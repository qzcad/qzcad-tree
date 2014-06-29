/**
  * @author Сергей Чопоров
  * @date 15/10/2012
  * @version 1.0.3
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
     * @brief Конструктор
     * @param x Начальное значение ординаты
     */
    Point1D(Floating x = 0);
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
     * @brief Ордината
     * @return Ординату точки или ее проекции в 1D
     */
    virtual Floating x() const;
    /**
     * @brief Абсцисса
     * @return 0.0
     */
    virtual Floating y() const;
    /**
     * @brief Аппликата
     * @return 0.0
     */
    virtual Floating z() const;
    /**
     * @brief Длина
     * @return Длину вектора (отрезка)
     */
    virtual Floating length() const;
    /**
     * @brief Установить координату точки
     * @param x Новое значение ординаты точки
     */
    void set(Floating x);
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
    Floating operator *(const Point1D &point) const;
    /// Произведение числа d на вектор, заданный точкой P
    friend Point1D operator *(Floating d, const Point1D &point);
    /// Отношение вектора, занного точкой point, к числу dec
    friend Point1D operator /(const Point1D &point, Floating dec);
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
    Floating distanceTo(const Point1D &point) const;
    /**
     * @brief Проверка на приблизительное равенство
     * @param point Точка, с которой производится сравнение
     * @param epsilon Точность
     * @return true, если точки приблизительно равны; false - иначе
     */
    bool isEqualTo(const Point1D &point, Floating epsilon = 0.00001) const;
private:
    Floating x_; //!< Ордината точки
};
}

#endif // POINT1D_H
