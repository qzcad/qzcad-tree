/**
  * @author Сергей Чопоров
  * @date 21/01/2014
  * @version 1.0.1
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef SEGMENTMESH2D_H
#define SEGMENTMESH2D_H

#include <vector>
#include <functional>
#include <list>

#include "mesh2d.h"
#include "point2d.h"
#include "segment.h"

namespace msh
{
/**
 * @brief Класс Path2D - контур в двумерном пространстве
 */
class SegmentMesh2D : public Mesh2D
{
public:
    /**
     * @brief Path2D Конструктор по умолчанию
     */
    SegmentMesh2D();
    /**
     * @brief Path2D Конструктор копирования
     * @param path2d Объект для копирования
     */
    SegmentMesh2D(const SegmentMesh2D &mesh);
    /**
     * @brief Конструктор создает копию объекта, переданого по ссылке
     * @param mesh Ссылка на объект для копирования
     */
    SegmentMesh2D(const SegmentMesh2D *mesh);
    /**
     * @brief Конструктор создает равномерную секту области, определенной функционально
     * @param xCount Количество узлов вдоль оси абсцисс
     * @param yCount Количество узлов вдоль оси ординат
     * @param xMin Абсцисса нижнего левого угла прямоугольной области
     * @param yMin Ордината нижнего левого угла прямоугольной области
     * @param width Ширина прямоугольной области
     * @param height Высота прямоугольной области
     * @param func Функция области
     * @param charPoint Список характерных точек
     */
    SegmentMesh2D(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double(double, double)> func, std::list<Point2D> charPoint);
    /**
     * @brief elementsCount Количество элементов
     * @return Количество отрезков (граней), которые образуют контур
     */
    UInteger elementsCount() const;
    /**
     * @brief Элемет сетки
     * @param number Номер элемента сетки
     * @return Указатель на элемент сетки
     */
    virtual ElementPointer element(const UInteger &number) const;
    /**
     * @brief Добавить элемент сетки
     * @param node0 Первый узел элемента
     * @param node1 Второй узел элемента
     */
    void addElement(const UInteger &node0, const UInteger &node1);
    /**
     * @brief Добавить элемент к сетке
     * @param Массив ссылок (номеров) на узлы
     */
    void addElement(const std::vector<UInteger> &nodes_ref);
    /**
     * @brief Изменить направление обхода во всех элементах сетки
     */
    virtual void directionChange();
    /**
     * @brief area Подсчитать площадь конечного элемента. Для данного типа элемента вычислется его длина
     * @param number Номер элемента
     * @return Площадь элемента (длина)
     */
    virtual double area(const UInteger &number) const;
    /**
     * @brief Определить принадлежность элемента границе
     * @param number Номер элемента
     * @return true - граничный элемент; false - внутренний
     */
    virtual bool isBorderElement(const UInteger &number) const;
    /**
     * @brief Отрезок с заданным номером
     * @param number Номер элемента (отрезка)
     * @return Отрезок с заданным номером
     */
    Segment segment(const UInteger &number) const;
protected:

private:
    std::vector<Segment> element_; //!< Массив элементов
};
}

#endif // SEGMENTMESH2D_H
