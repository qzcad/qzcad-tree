/**
  * @author Сергей Чопоров
  * @date 12/02/2014
  * @version 1.0.2
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef QUADRILATERALMESH2D_H
#define QUADRILATERALMESH2D_H

#include <functional>
#include <list>

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
     * @brief Метод создает равномерную структурированную секту в прямоугольной области
     * @param xCount Количество узлов вдоль оси абсцисс
     * @param yCount Количество узлов вдоль оси ординат
     * @param xMin Абсцисса нижнего левого угла прямоугольной области
     * @param yMin Ордината нижнего левого угла прямоугольной области
     * @param width Ширина прямоугольной области
     * @param height Высота прямоугольной области
     */
    void rectangleDomain(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height);
    /**
     * @brief Метод создает равномерную структурированную сетку для выпуклой четырехугольной области
     * @param xCount Количество узлов по первому направлению (ребра 0-1 и 2-3)
     * @param yCount Количество узлов по второму направлению (ребра 1-2 и 3-0)
     * @param v0 Координаты узла 0
     * @param v1 Координаты узла 1
     * @param v2 Координаты узла 2
     * @param v3 Координаты узла 3
     */
    void quadDomain(const UInteger &xCount, const UInteger &yCount, const Point2D &v0, const Point2D &v1, const Point2D &v2, const Point2D &v3);
    /**
     * @brief Метод создает сетку для треугольной области
     * @param count Количество узлов на сторону треугольника (должно быть четным)
     * @param v0 Координаты узла 0
     * @param v1 Координаты узла 1
     * @param v2 Координаты узла 2
     */
    void triangleDomain(const UInteger &count, const Point2D &v0, const Point2D &v1, const Point2D &v2);
    /**
     * @brief Метод создает блочно-структурированную сетку для круга (части круга)
     * @param count Базовое количество узлов (По окружности будет 4n - для целого круга, 4n - для половинки, 2n - для четверти)
     * @param center Координаты центра
     * @param radius Радиус
     * @param part Часть круга для дискретизации (возможные значение: 1 - целый круг, 2 - половинка, 4 - четверть)
     */
    void circleDomain(const UInteger &count, const Point2D &center, const double &radius, unsigned short part = 1);
    /**
     * @brief Конструктор копирования
     * @param mesh Экземпляр объекта для копирования
     */
    QuadrilateralMesh2D(const QuadrilateralMesh2D &mesh);
    /**
     * @brief Конструктор создает копию объекта, переданого по ссылке
     * @param mesh Ссылка на объект для копирования
     */
    QuadrilateralMesh2D(const QuadrilateralMesh2D *mesh);
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
    void functionalDomain(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double(double, double)> func, std::list<Point2D> charPoint);
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
    virtual double area(const UInteger &number) const;
    /**
     * @brief Добавить элемент к  сетке
     * @param node0 Номер узла в вершине 0
     * @param node1 Номер узла в вершине 1
     * @param node2 Номер узла в вершине 2
     * @param node3 Номер узла в вершине 3
     * Этот метод не является потокобезопасным.
     */
    void addElement(const UInteger &node0, const UInteger &node1, const UInteger &node2, const UInteger &node3);
    /**
     * @brief Добавить элемент к сетке
     * @param quad Ссылка на элемент, который необходимо добавить к сетке
     * Этот метод не является потокобезопасным.
     */
    void addElement(const Quadrilateral &quad);
    /**
     * @brief Добавить элемент к сетке
     * @param Массив ссылок (номеров) на узлы
     */
    void addElement(const std::vector<UInteger> &nodes_ref);
    /**
     * @brief Метод очищает информацию об елементах
     */
    virtual void clearElements();
protected:
    /**
     * @brief Функция формы изопараметрического четырехугольного элемента
     * @param i Номер узла
     * @param xi Значение параметра первого направления
     * @param eta Значение параметра второго направления
     * @return Значение функции формы
     */
    double isoFunc(const UInteger &i, const double &xi, const double &eta);

    double functional(const std::vector<double> &vars);
    /**
     * @brief Метод для добавления сеток на основе трансфинитной интерполяции
     * Все узлы добавляются как внутренние (INNER, @see NodeType). Функции top, bottom, left, right - функторы (functors), указатели на функции (function pointers) или лямбды (lambda), которые возвращают координаты точки кривой в соответствии с заданным параметром.
     * Ограничения: top(0) = left(1), top(1) = right(1), bottom(0) = left(0), bottom(1) = right(0)
     * @param top Верхняя граница (функция xi)
     * @param bottom Нижняя граница (функция xi)
     * @param left Левая граница (функция eta)
     * @param roght Правая граница (функция eta)
     * @param xiCount Количество узлов по "направлению" xi
     * @param etaCount Количество узлов по "направлению" eta
     */
    template<typename TopFunc, typename BottomFunc, typename LeftFunc, typename RightFunc>
    void addTransfiniteMesh(TopFunc top, BottomFunc bottom, LeftFunc left, RightFunc right, const UInteger &xiCount, const UInteger &etaCount);
protected:
    std::vector<Quadrilateral> element_; //!< Массив элементов
};
}

#endif // QUADRILATERALMESH2D_H
