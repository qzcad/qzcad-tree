/**
  * @author Сергей Чопоров
  * @date 18/06/2014
  * @version 1.0.1
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef POINT3D_H
#define POINT3D_H

#include "point2d.h"

namespace msh {
/**
 * @brief Класс Point3D - точка в пространстве (3D)
 * @see Point2D
 */
class Point3D : public Point2D
{
public:
    /**
     * @brief Конструктор по умолчанию инициализирует точку началом координат (0; 0; 0)
     */
    Point3D();
    /**
     * @brief Конструктор
     * @param x Абсцисса
     * @param y Ордината
     * @param z Аппликата
     */
    Point3D(const double &x, const double &y, const double &z);
    /**
     * @brief Конструктор копирования
     * @param point
     */
    Point3D(const Point3D &point);
    /**
     * @brief Конструктор, который формирует точку как вершину вектора, задоного двумя точками (разность между первой и второй точками)
     * @param firstPoint Первая точка
     * @param secondPoint Вторая точка
     */
    Point3D(const Point3D &firstPoint, const Point3D &secondPoint);
    /**
     * @brief Размерность пространства
     * @return размерность пространства, в ктором определена точка
     */
    virtual int dimension() const;
    /**
     * @brief Аппликата
     * @return Аппликату точки или ее проекции
     */
    virtual double z() const;
    /**
     * @brief Длина вектора
     * @return Длину вектора, определенного началом координат и координатами точки
     */
    virtual double length() const;
    /**
     * @brief Координаты нормализованного вектора
     * @return Координаты нормализованного вектора (единичного вектора того же направления)
     */
    Point3D normalized() const;
    /**
     * @brief Метод для вычисления расстояния до точки
     * @param point Точка, до которой необходимо определить расстояние
     * @return Расстояние до точки point
     */
    double distanceTo(const Point3D &point) const;
    /**
     * @brief Метод для вычисления расстояния от точки до отрезка
     * @param segment0 Начало отрезка
     * @param segment1 Конец отрезка
     * @return Расстояние до ближайшей точки отрезка
     */
    double distanceTo(const Point3D &segment0, const Point3D &segment1) const;
    /**
     * @brief Метод для вычисления расстояние от точки до треугольника
     * Вершины треугольника должны быть переданы в порядке обхода против часовой стрелки
     * @param triangle0 Первая вершина треугольника
     * @param triangle1 Вторая вершина треугольника
     * @param triangle2 Третья вершина треугольника
     * @return Расстояние до ближайшей точки треугольника
     */
    double distanceTo(const Point3D &triangle0, const Point3D &triangle1, const Point3D &triangle2) const;
    /**
     * @brief Проверка на приблизительное равенство
     * @param point Точка, с которой производится сравнение
     * @param epsilon Точность
     * @return true, если точки приблизительно равны; false - иначе
     */
    bool isEqualTo(const Point3D &point, double epsilon = 0.00001) const;
    /**
     * @brief Векторное произведение векторов, определенных текущей точкой и заданной
     * @param point Координаты ветора, с которым необходимо найти векторное произведение
     * @return Значение векторного произведения
     */
    Point3D product(const Point3D &point) const;
    /**
     * @brief Установить координаты точки
     * @param x Абсцисса
     * @param y Ордината
     * @param z Аппликата
     */
    void set(const double &x, const double &y, const double &z);
    /**
     * @brief Установить значение аппликаты
     * @param z Новое значение аппликаты
     */
    void setZ(const double &z);
    /**
     * @brief Напечатать на стандартную консоль координаты
     */
    void print() const;
    /**
     * @brief Напечатать на стандартную консоль координаты с переходом на новую строку
     */
    void println() const;
    /**
     * @brief Получить координаты точки в новой ортогональной системе координат, определенной тройкой векторов
     * @param Vx Вектор первого направления
     * @param Vy Вектор второго напраления
     * @param Vz Вектор третьего напраления
     * @return Координаты точки в заданной системе
     */
    Point3D inCoordSystem(const Point3D &Vx, const Point3D &Vy, const Point3D &Vz) const;
    /**
     * @brief Получить координаты точки, заданной в системе координат Vx, Vy, Vz, в декартовой системе
     * @param Vx Вектор первого направления
     * @param Vy Вектор второго напраления
     * @param Vz Вектор третьего напраления
     * @return Координаты точки в декартовой системе
     */
    Point3D fromCoordSystem(const Point3D &Vx, const Point3D &Vy, const Point3D &Vz) const;
    /**
     * @brief Масштабировать координаты точки на коэффициент
     * @param d Коэффициент масштабирования
     */
    virtual void scale(const double &d);
    /**
     * @brief Вычислить угол B-this-C
     * @param B Первая точка угла
     * @param C Вторая точка угла
     * @return Угол B-this-C в радианах
     */
    double angle(const Point3D &B, const Point3D &C) const;

    /// @name Операторы
    /// @{
    /// Присваивание
    Point2D &operator =(const Point3D &point);
    /// Проверка на равенство leftPoint и rightPont
    friend bool operator ==(const Point3D &leftPoint, const Point3D &rightPoint);
    /// Унарный минус
    friend const Point3D operator -(const Point3D &point);
    /// Бинарный минус векторов, определенных точками leftPoint и rightPoint
    friend const Point3D operator -(const Point3D &leftPoint, const Point3D &rightPoint);
    /// Бинарный плюс векторов, определенных точками leftPoint и rightPoint
    friend const Point3D operator +(const Point3D &leftPoint, const Point3D &rightPoint);
    /// Скалярное произведение векторов, определенных точками this и point
    double operator *(const Point3D &point) const;
    /// Произведение числа dec на вектор, заданный точкой point
    friend const Point3D operator *(double dec, const Point3D &point);
    /// Отношение вектора, занного точкой point, к числу dec
    friend const Point3D operator /(const Point3D &point, double dec);
    /**
     * @brief Оператор сравнения "меньше" с использованием лексикографического порядка
     * @param leftPoint Левая точка
     * @param rightPoint Правая точка
     * @return leftPoint < rightPoint
     */
    friend bool operator <(const Point3D &leftPoint, const Point3D &rightPoint);
    /**
     * @brief Оператор сравнения "меньше или равно" с использованием лексикографического порядка
     * @param leftPoint Левая точка
     * @param rightPoint Правая точка
     * @return leftPoint < rightPoint
     */
    friend bool operator <=(const Point3D &leftPoint, const Point3D &rightPoint);
    /// @}
    /**
     * @brief Проверка скрещиваются ли прямые P0P1 и Q0Q1
     * @param P0 Первая точка прямой P
     * @param P1 Вторая точка прямой P
     * @param Q0 Первая точка прямой Q
     * @param Q1 Вторая точка прямой Q
     * @param p Коэффицент, соответсвующий точке пресечения на прямой P ( X = P0 + p * (P1 - P0) )
     * @param q Коэффицент, соответсвующий точке пресечения на прямой q ( X = Q0 + q * (Q1 - Q0) )
     * @return
     */
    friend bool isSkew(const Point3D &P0, const Point3D &P1, const Point3D &Q0, const Point3D &Q1, double &p, double &q);
    /**
     * @brief Вычислить внешнюю нормаль к треугольнику ABC
     * @param A Первая вершина треугольника
     * @param B Вторая вершина треугольника
     * @param C Третья вершина треугольника
     * @return Нормализовнные координаты вектора-нормали
     */
    friend Point3D normal3(const Point3D &A, const Point3D &B, const Point3D &C);
private:
    double z_; //!< Аппликата точки
};
}

#endif // POINT3D_H
