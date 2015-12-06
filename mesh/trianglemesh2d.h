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
#include <list>

#include "mesh2d.h"
#include "triangle.h"
#include "segmentmesh2d.h"

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
     * @param charPoint Список характерных точек
     */
    TriangleMesh2D(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double(double, double)> func, std::list<Point2D> charPoint);
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
     * @brief Конструктор для построения триангуляции Делоне заданного двумерного конутра
     * @param mesh Указатель на контур
     */
    TriangleMesh2D(const SegmentMesh2D *mesh);
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
    /**
     * @brief Добавить элемент к сетке
     * @param Массив ссылок (номеров) на узлы
     */
    void addElement(const std::vector<UInteger> &nodes_ref);
    /**
     * @brief Вычислить значение якобиана элемента (линейный случай)
     * @param elementNum Номер элемента
     * Для линейных функций формы треугольного элемента значение якобиана не зависит от параметров xi (L1) и eta (L2).
     * @return Значение якобиана элемента номер elementNum
     */
    double jacobian(const UInteger &elementNum);
    /**
     * @brief Вычислить соотношение длин сторон (минимальной к максимальной)
     * @param elNum Номер элемента
     * @return Соотношение длин сторон (минимальной к максимальной)
     */
    double lengthAspect(const UInteger &elNum);
    /**
     * @brief Вычислить значение минимального угла в элементе
     * @param elNum Номер элемента
     * @return Минимальный угол элемента (радианы)
     */
    double minAngle(const UInteger &elNum);
    /**
     * @brief Вычислить соотношение углов элемента
     * @param elNum Номер элемента
     * @return Соотношение углов элемента (радианы)
     */
    double angleAspect(const UInteger &elNum);
    /**
     * @brief Триангуляция Делоне объекта, представленного функционально
     * @param xCount Количество узлов вдоль оси абсцисс
     * @param yCount Количество узлов вдоль оси ординат
     * @param xMin Абсцисса нижнего левого угла прямоугольной области
     * @param yMin Ордината нижнего левого угла прямоугольной области
     * @param width Ширина прямоугольной области
     * @param height Высота прямоугольной области
     * @param func Функция области
     * @param charPoint Список характерных точек
     */
    void delaunay(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double(double, double)> func, std::list<Point2D> charPoint);
    /**
     * @brief Триангуляция Делоне с использованием метода Рапперта для объекта, представленного функционально
     * @param xCount Количество узлов вдоль оси абсцисс
     * @param yCount Количество узлов вдоль оси ординат
     * @param xMin Абсцисса нижнего левого угла прямоугольной области
     * @param yMin Ордината нижнего левого угла прямоугольной области
     * @param width Ширина прямоугольной области
     * @param height Высота прямоугольной области
     * @param func Функция области
     * @param charPoint Список характерных точек
     * @param refineArea Флаг, указывающий на необходимость оптимизации по площади
     */
    void ruppert(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double(double, double)> func, std::list<Point2D> charPoint, bool refineArea = false);
    /**
     * @brief Элемент сетки (треугольник)
     * @param number Номер треугольника
     * @return Треугольник с указанным номером
     */
    Triangle triangle(const UInteger &number) const;
protected:
    /**
     * @brief Метод находит значение минимального угла в треугольнике, определенном координатами вершин
     * @param A Координаты вершины
     * @param B Координаты вершины
     * @param C Координаты вершины
     *       C
     *       ^
     *     b/ \a
     *    A --- B
     *       c
     * @return Значение минимального угла в треугольнике (радианы)
     */
    double minAngle(const Point2D &A, const Point2D &B, const Point2D &C);
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
    bool angles(const Point2D &A, const Point2D &B, const Point2D &C, double &alpha, double &beta, double &gamma);
    /**
     * @brief Структура для временного хранения триангуляции
     */
    struct Triangulation
    {
        std::vector<Point2D> nodes;
        std::vector<NodeType> types;
        std::list<Triangle> triangles;
    };
    /**
     * @brief Метод для построения триангуляции Делоне для заданного контура, используя суперобласть
     * @param mesh Указатель на сетку-контур
     * @return Триаунгуляцию, как выпуклого множества
     */
    Triangulation superDelaunay(const SegmentMesh2D *mesh);
    /**
     * @brief Оптимизация триангуляции Делоне методом Рапперта (супер область)
     * @param Triangulation Ссылка на триангуляцию Делоне супер области
     */
    void superRuppert(Triangulation &triangulation);
    /**
     * @brief Метод контроля за максимальной площадью элемента на основе вставки нового узла в центр масс
     * @param func Функция области (обрабатываются только внутренние элементы)
     * @param triangulation Ссылка на триангуляцию Делоне супер области
     */
    void areaRefinement(double max_area, std::function<double(double, double)> func, Triangulation &triangulation);
    /**
     * @brief Операция вставки узла в триангуляцию Делоне
     * @param point Координаты узла для вставки (входной парметр)
     * @param nodes Ссылка на массив узлов для вставки (выходной параметр)
     * @paran triangles Ссылка на массив элементов для вставки (выходной параметр)
     */
    bool insertDelaunayNode(const Point2D &point, const NodeType &type, std::vector<Point2D> &nodes, std::vector<NodeType> &types, std::list<Triangle> &triangles);
    /**
     * @brief Метод вычисления площади треугольника с учетом знака (площадь отрицательная при обходе узлов треугольника против часовой стрелки)
     * @param A Координаты первого узла треугольника
     * @param B Координаты второго узла треугольника
     * @param C Координаты третьего узла треугольника
     * @return Площадь треугольника со знаком
     */
    double signedArea(const Point2D &A, const Point2D &B, const Point2D &C) const;
    /**
     * @brief Метод проверки вхождения точки в описанную вокруг треугольника окргужность
     * @param xp Абсцисса точки
     * @param yp Ордината точки
     * @param x1 Абсцисса первого узла треугольника
     * @param y1 Ордината первого узла треугольника
     * @param x2 Абсцисса второго узла треугольника
     * @param y2 Ордината второго узла треугольника
     * @param x3 Абсцисса третьего узла треугольника
     * @param y3 Ордината третьего узла треугольника
     * @param xc Абсцисса центра описанной окружности
     * @param yc Ордината центра описанной окружности
     * @param r Радиус описанной окружности
     * @return true, если заданная точка попала в описанную окружность
     */
    bool circumCircle(double xp, double yp, double x1, double y1, double x2, double y2, double x3, double y3, double &xc, double &yc, double &r);
protected:
    std::vector<Triangle> element_; //!< Массив элементов
};
}

#endif // TRIANGLEMESH2D_H
