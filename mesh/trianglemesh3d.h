/**
  * @author Сергей Чопоров
  * @date 06/08/2015
  * @version 1.0.0
  * @copyright Copyright 2015 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef TRIANGLEMESH3D_H
#define TRIANGLEMESH3D_H

#include <functional>
#include <list>

#include "mesh3d.h"
#include "triangle.h"

namespace msh {
/**
 * @brief Класс TriangleMesh3D - дискртеная модель оболчки на базе треугольников
 */
class TriangleMesh3D : public Mesh3D
{
public:
    /**
     * @brief Конструктор по умолчанию
     */
    TriangleMesh3D();
    /**
     * @brief Конструктор копирования
     * @param mesh Экземпляр объекта для копирования
     */
    TriangleMesh3D(const TriangleMesh3D &mesh);
    /**
     * @brief Конструктор создает копию объекта, переданого по ссылке
     * @param mesh Ссылка на объект для копирования
     */
    TriangleMesh3D(const TriangleMesh3D *mesh);
    /**
     * @brief Конструктор для создания равномерной сетки в цилиндрических координатах
     * @param rCount Количество элементов вдоль радиуса
     * @param lCount Количество элементов вдоль образующей
     * @param radius Радиус
     * @param length Высота
     */
    TriangleMesh3D(const UInteger &rCount, const UInteger &lCount, const double &radius, const double &length);
    /**
     * @brief Конструктор для создания равномерной сетки в цилиндрических координатах с адаптацией под область, ограниченную функцией
     * @param rCount Количество элементов вдоль радиуса
     * @param lCount Количество элементов вдоль образующей
     * @param radius Радиус
     * @param length Высота
     * @param func Функция, положительная во внутренних точка и отрицательная во внешних
     * @param charPoints Массив характерных точек
     */
    TriangleMesh3D(const UInteger &rCount, const UInteger &lCount, const double &radius, const double &length, std::function<double(double, double, double)> func, std::list<Point3D> charPoints);
    /**
     * @brief Конструктор для создания равномерной сетки в конических координатах
     * @param rCount Количество элементов вдоль радиуса
     * @param lCount Количество элементов вдоль образующей
     * @param bottom_radius Нижний радиус
     * @param top_radius Верхний радиус
     * @param length Высота
     */
    TriangleMesh3D(const UInteger &rCount, const UInteger &lCount, const double &bottom_radius, const double &top_radius, const double &length);
    /**
     * @brief Количество элементов
     * @return Количество элементов в сетке
     */
    virtual UInteger elementsCount() const;
    /**
     * @brief Элемет сетки
     * @param number Номер элемента сетки
     * @return Указатель на копию элемента сетки
     */
    virtual ElementPointer element(const UInteger &number) const;
    /**
     * @brief Определить принадлежность элемента границе
     * @param number Номер элемента
     * @return true - граничный элемент; false - внутренний
     */
    virtual bool isBorderElement(const UInteger &number) const;
    /**
     * @brief Вычислить площадь поверхности дискретной модели
     * @return Площадь поверхности дискретной модели
     */
    virtual double surfaceArea() const;
    /**
     * @brief Добавить элемент к  сетке
     * @param node0 Номер узла в вершине 0
     * @param node1 Номер узла в вершине 1
     * @param node2 Номер узла в вершине 2
     * Этот метод не является потокобезопасным.
     */
    void addElement(const UInteger &node0, const UInteger &node1, const UInteger &node2);
    /**
     * @brief Добавить элемент к сетке
     * @param triangle Ссылка на элемент, который необходимо добавить к сетке
     * Этот метод не является потокобезопасным.
     */
    void addElement(const Triangle &triangle);
    /**
     * @brief area Подсчитать площадь элемента
     * @param number Номер элемента
     * @return Площадь элемента
     */
    virtual double area(const UInteger &number) const;
    /**
     * @brief Вычислить значение минимального угла в элементе
     * @param elNum Номер элемента
     * @return Минимальный угол элемента (радианы)
     */
    double minAngle(const UInteger &elNum);
protected:
    /**
     * @brief Функция для подсчета значений углов треугольника
     * @param A Координаты первой вершины
     * @param B Координаты второй вершины
     * @param C Координаты третей вершины
     * @param alpha Угол в вершине A (радианы)
     * @param beta Угол в вершине B (радианы)
     * @param gamma Угол в вершине C (радианы)
     * @return true, если трейгольник невырожденный, иначе - false
     *       C
     *       ^
     *     b/ \a
     *    A --- B
     *       c
     */
    bool angles(const Point3D &A, const Point3D &B, const Point3D &C, double &alpha, double &beta, double &gamma);
    /**
     * @brief Метод находит значение минимального угла в треугольнике, определенном координатами вершин
     * @param A Координаты вершины
     * @param B Координаты вершины
     * @param C Координаты вершины
     * @return Значение минимального угла в треугольнике (радианы)
     */
    double minAngle(const Point3D &A, const Point3D &B, const Point3D &C);
protected:
    std::vector<Triangle> element_; //!< Массив элементов
};
}
#endif // TRIANGLEMESH3D_H
