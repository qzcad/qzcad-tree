/**
  * @author Сергей Чопоров
  * @date 13/04/2015
  * @version 1.0.0
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef TRIANGLEMESH2D_H
#define TRIANGLEMESH2D_H

#include <functional>

#include "mesh2d.h"
#include "triangle.h"

namespace msh
{
/**
 * @brief Класс TriangleMesh2D - абстракция сетки треугольных элементов на плоскости.
 */
class TriangleMesh2D : public Mesh2D
{
public:
    /**
     * @brief Коснтруктор по умолчанию.
     * Создается пустая сетка.
     */
    TriangleMesh2D();
    /**
     * @brief Конструктор создает равномерную структурированную секту в прямоугольной области
     * @param xCount Количество узлов вдоль оси абсцисс
     * @param yCount Количество узлов вдоль оси ординат
     * @param xMin Абсцисса нижнего левого угла прямоугольной области
     * @param yMin Ордината нижнего левого угла прямоугольной области
     * @param width Ширина прямоугольной области
     * @param height Высота прямоугольной области
     */
    TriangleMesh2D(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height);
    /**
     * @brief Конструктор создает равномерную секту области, определенной функционально
     * @param xCount Количество узлов вдоль оси абсцисс
     * @param yCount Количество узлов вдоль оси ординат
     * @param xMin Абсцисса нижнего левого угла прямоугольной области
     * @param yMin Ордината нижнего левого угла прямоугольной области
     * @param width Ширина прямоугольной области
     * @param height Высота прямоугольной области
     * @param func Функция области
     */
    TriangleMesh2D(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double(double, double)> func);
    /**
     * @brief Конструктор копирования
     * @param mesh Экземпляр объекта для копирования
     */
    TriangleMesh2D(const TriangleMesh2D &mesh);
    /**
     * @brief Конструктор создает копию объекта, переданого по ссылке
     * @param mesh Ссылка на объект для копирования
     */
    TriangleMesh2D(const TriangleMesh2D *mesh);
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
     * @brief Изменить направление обхода во всех элементах сетки
     */
    virtual void directionChange();
    /**
     * @brief area Подсчитать площадь конечного элемента
     * @param number Номер элемента
     * @return Площадь элемента
     */
    virtual double area(const UInteger &number) const;
    /**
     * @brief Добавить элемент (треугольник)
     * @param node0 Номер (хеш-код) узла
     * @param node1 Номер (хеш-код) узла
     * @param node2 Номер (хеш-код) узла
     * Этот метод не является потокобезопасным.
     */
    void addElement(const UInteger &node0, const UInteger &node1, const UInteger &node2);
protected:
    std::vector<Triangle> element_; //!< Массив элементов
};
}

#endif // TRIANGLEMESH2D_H
