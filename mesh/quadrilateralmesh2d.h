/**
  * @author Сергей Чопоров
  * @date 12/02/2014
  * @version 1.0.2
  * 1.0.2:
  *     - добавлен конструктор по умолчанию
  */
#ifndef QUADRILATERALMESH2D_H
#define QUADRILATERALMESH2D_H

#include "mesh2d.h"
#include "quadrilateral.h"

namespace msh
{
/**
 * @brief Класс QuadrilateralMesh2D - неструктурированная сетка четырехугольных элементов
 */
class QuadrilateralMesh2D : public Mesh2D
{
public:
    /**
     * @brief Коснтруктор по умолчанию
     */
    QuadrilateralMesh2D();
    /**
     * @brief Конструктор создает равномерную структурированную секту в прямоугольной областе
     * @param xCount Количество узлов вдоль оси ординат
     * @param yCount Количество узлов вдоль оси абсцисс
     * @param xMin Ордината нижнего левого угла прямоугольной области
     * @param yMin Абсцисса нижнего левого угла прямоугольной области
     * @param width Ширина прямоугольной области
     * @param height Высота прямоугольной области
     */
    QuadrilateralMesh2D(const UInteger &xCount, const UInteger &yCount, const Floating &xMin, const Floating &yMin, const Floating &width, const Floating &height);
    /**
     * @brief Конструктор создает равномерную структурированную сетку для выпуклой четырехугольной области
     * @param xCount Количество узлов по первому направлению (ребра 0-1 и 2-3)
     * @param yCount Количество узлов по второму направлению (ребра 1-2 и 3-0)
     * @param v0 Координаты узла 0
     * @param v1 Координаты узла 1
     * @param v2 Координаты узла 2
     * @param v3 Координаты узла 3
     */
    QuadrilateralMesh2D(const UInteger &xCount, const UInteger &yCount, const Point2D &v0, const Point2D &v1, const Point2D &v2, const Point2D &v3);
    /**
     * @brief Конструктор создает сетку для треугольной области
     * @param count Количество узлов на сторону треугольника (должно быть четным)
     * @param v0 Координаты узла 0
     * @param v1 Координаты узла 1
     * @param v2 Координаты узла 2
     */
    QuadrilateralMesh2D(const UInteger &count, const Point2D &v0, const Point2D &v1, const Point2D &v2);
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
    void minimizeFunctional();
    virtual Floating nodeValue(const UInteger &number) const;
    /**
     * @brief Определить принадлежность элемента границе
     * @param number Номер элемента
     * @return true - граничный элемент; false - внутренний
     */
    virtual bool isBorderElement(const UInteger &number) const;
    /**
     * @brief Четырехугольник, соответсвующий заданному номеру
     * @param number Номер четырехугольника
     * @return Четырехугольник с номером number
     */
    Quadrilateral quadrilateral(const UInteger &number) const;
    /**
     * @brief Изменить направление обхода во всех элементах сетки
     */
    virtual void directionChange();
    /**
     * @brief area Подсчитать площадь элемента
     * @param number Номер элемента
     * @return Площадь элемента
     */
    virtual Floating area(const UInteger &number);
protected:
    /**
     * @brief Добавить элемент к  сетке
     * @param node0 Номер узла в вершине 0
     * @param node1 Номер узла в вершине 1
     * @param node2 Номер узла в вершине 2
     * @param node3 Номер узла в вершине 3
     */
    void addElement(const UInteger &node0, const UInteger &node1, const UInteger &node2, const UInteger &node3);
    /**
     * @brief Функция формы изопараметрического четырехугольного элемента
     * @param i Номер узла
     * @param xi Значение параметра первого направления
     * @param eta Значение параметра второго направления
     * @return Значение функции формы
     */
    Floating isoFunc(const UInteger &i, const Floating &xi, const Floating &eta);

    Floating functional(Floating *vars, const std::vector<UInteger> &nodeVariable);
    void nabla(const UInteger &size, Floating *x,
               const std::vector<UInteger> &nodeVariable, Floating *gradient, Floating h = 0.0001);
    Floating lambda(const UInteger &size, Floating *x, Floating *s, const Floating &lambda_val,
                    const std::vector<UInteger> &nodeVariable);
    Floating goldenRatio(const UInteger &size, const Floating &a, const Floating &b, Floating *x0, Floating *s, const std::vector<UInteger> &nodeVariable, Floating epsilon = 0.0001, UInteger maxIter = 300);
    Floating norm2(const UInteger &size, Floating *x);
    void conjugateGradient(const UInteger &size, Floating *x0, Floating *xMin,
                           const std::vector<UInteger> &nodeVariable, Floating epsilon = 0.0001, UInteger maxIter = 300);
//    double sqr(const Floating &a) { return a * a; }
protected:
    std::vector<Quadrilateral> element_; //!< Массив элементов
};
}

#endif // QUADRILATERALMESH2D_H
