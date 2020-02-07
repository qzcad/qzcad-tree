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
#include "quadrilateralmesh2d.h"

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
    void fromQuadrilateralMesh(const QuadrilateralMesh2D *mesh);
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
     * @brief Метод создает равномерную секту области, определенной функционально
     * @param xCount Количество узлов вдоль оси абсцисс
     * @param yCount Количество узлов вдоль оси ординат
     * @param xMin Абсцисса нижнего левого угла прямоугольной области
     * @param yMin Ордината нижнего левого угла прямоугольной области
     * @param width Ширина прямоугольной области
     * @param height Высота прямоугольной области
     * @param func Функция области
     * @param charPoint Список характерных точек
     * @param smooth Количество итераций сглаживания
     * @param optimize Количество итераций оптимизации
     */
    void functionalDomain(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double(double, double)> func, std::list<Point2D> charPoint, int smooth = 0, int optimize = 0);
    void backgroundGrid(const TriangleMesh2D *mesh, std::function<double(double, double)> primary, std::list<Point2D> charPoint, double level = 0.0, int smooth = 0, int optimize = 0);
    /**
     * @brief Метод построения модели с использованием фоновых треугольников для области, определеной парой функций
     * @param primary Указатель на функцию, которая определяет границы области
     * @param secondary Указатель на функцию для адаптации узлов (например, геометрия слоя материала)
     * @param charPoint Список характерных точек
     * @param level Значение функции, соответствующее границе
     * @param smooth Количество итераций сглаживания
     * @param optimize Количество итераций оптимизации
     */
    void backgroundGrid(const TriangleMesh2D *mesh, std::function<double(double, double)> primary, std::function<double(double, double)> secondary, std::list<Point2D> charPoint, double level = 0.0, int smooth = 0, int optimize = 0);
    /**
     * @brief Конструктор копирования
     * @param mesh Экземпляр объекта для копирования
     */
    TriangleMesh2D(const TriangleMesh2D &mesh);
    /**
     * @brief Конструктор создает копию объекта, переданого по указателю
     * @param mesh Указатель на объект для копирования
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
     * @brief Добавить элемент к сетке
     * @param triangle Экземпляр элемента (треугольник)
     */
    void addElement(const Triangle &triangle);
    /**
     * @brief Добавить элемент к сетке
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
    double jacobian(const UInteger &elementNum) const;
    /**
     * @brief Вычислить значение минимального угла в элементе
     * @param elNum Номер элемента
     * @return Минимальный угол элемента (радианы)
     */
    virtual double minAngle(const UInteger &elNum) const;
    /**
     * @brief Вычислить значение максимального угла в элементе
     * @param elNum Номер элемента
     * @return Минимальный угол элемента (радианы)
     */
    virtual double maxAngle(const UInteger &elNum) const;
    /**
     * @brief Вычислить соотношение углов элемента
     * @param elNum Номер элемента
     * @return Соотношение углов элемента (радианы)
     */
    double angleAspect(const UInteger &elNum) const;
    /**
     * @brief Триангуляция Делоне контура функциональной модели
     * @param mesh Объект-конутр
     * @param func Указатель на функцию области. Если он равен nullptr, то используется сfunction.
     * @see SegmentMesh2D::cfunction
     */
    void delaunay(const SegmentMesh2D &mesh, std::function<double(double, double)> func);
    /**
     * @brief Триангуляция Делоне с использованием улучшения модели методом Рапперта
     * @param mesh Объект-контур
     * @param func Указатель на функцию области. Если он равен nullptr, то используется сfunction.
     * @param alpha Минимально допустимый угол в радианах
     * @param max_area Если больше нуля, то удаляются все треугольники, площадь которых больше заданного значения.
     */
    void ruppert(const SegmentMesh2D &mesh, std::function<double(double, double)> func, double alpha = 0.436332, double max_area = 0.0);
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
     * @brief Триангуляция Делоне с использованием метода Рапперта для контактных задач
     * @param xCount Количество узлов вдоль оси абсцисс
     * @param yCount Количество узлов вдоль оси ординат
     * @param xMin Абсцисса нижнего левого угла прямоугольной области
     * @param yMin Ордината нижнего левого угла прямоугольной области
     * @param width Ширина прямоугольной области
     * @param height Высота прямоугольной области
     * @param func_a Функция области A
     * @param func_b Функция области B
     * @param charPoint Список характерных точек
     * @param delta Параметр сгущения элементов в окрестности контакта (если меньше 0, то сгущение отсутствует)
     * @param refineArea Флаг, указывающий на необходимость оптимизации по площади
     */
    void ruppert(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double(double, double)> func_a, std::function<double(double, double)> func_b, std::list<Point2D> charPoint, double delta = -1.0, bool refineArea = false);
    /**
     * @brief Элемент сетки (треугольник)
     * @param number Номер треугольника
     * @return Треугольник с указанным номером
     */
    Triangle triangle(const UInteger &number) const;
    /**
     * @brief Метод очищает информацию об елементах
     */
    virtual void clearElements();
    /**
     * @brief Метод реализует операцию смены диагонали между парой треугольников
     * @param print_messages если true, то будет выводится на экран информация о ходе выполнения
     */
    void flip(bool print_messages=true);
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
    static double minAngle(const Point2D &A, const Point2D &B, const Point2D &C);
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
    static bool angles(const Point2D &A, const Point2D &B, const Point2D &C, double &alpha, double &beta, double &gamma);
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
    static Triangulation superDelaunay(SegmentMesh2D *mesh, std::function<double(double, double)> func);
    /**
     * @brief Оптимизация триангуляции Делоне методом Рапперта (супер область)
     * @param Triangulation Ссылка на триангуляцию Делоне супер области
     */
    static void superRuppert(Triangulation &triangulation, SegmentMesh2D *mesh, std::function<double(double, double)> func, double alpha = 0.436332);
    /**
     * @brief Метод выполняет процедуру разбиения всех сегментов, в диаметр-окружности которых поападают другие вершины
     * @param triangulation Ссылка на триангуляцию Делоне
     */
    static void splitSegments(Triangulation &triangulation);
    /**
     * @brief Метод контроля за максимальной площадью элемента на основе вставки нового узла в центр масс
     * @param func Функция области (обрабатываются только внутренние элементы)
     * @param triangulation Ссылка на триангуляцию Делоне супер области
     */
    static void areaRefinement(double max_area, std::function<double(double, double)> func, Triangulation &triangulation);
    /**
     * @brief Операция вставки узла в триангуляцию Делоне
     * @param point Координаты узла для вставки (входной парметр)
     * @param type Тип добавляемого узла (входной парметр)
     * @paran triangles Ссылка на массив элементов для вставки (выходной параметр)
     */
    static bool insertDelaunayNode(const Point2D &point, const NodeType &type, Triangulation &triangulation);
    /**
     * @brief Метод вычисления площади треугольника с учетом знака (площадь отрицательная при обходе узлов треугольника против часовой стрелки)
     * @param A Координаты первого узла треугольника
     * @param B Координаты второго узла треугольника
     * @param C Координаты третьего узла треугольника
     * @return Площадь треугольника со знаком
     */
    static double signedArea(const Point2D &A, const Point2D &B, const Point2D &C);
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
    static bool circumCircle(const double &xp,
                             const double &yp,
                             const double &x1,
                             const double &y1,
                             const double &x2,
                             const double &y2,
                             const double &x3,
                             const double &y3,
                             double &xc,
                             double &yc,
                             double &r);
    /**
     * @brief Разбить элементы (сгустить сетку)
     * @param eNumbers Список номеров элементов
     */
    void subdivide(std::list<UInteger> eNumbers, std::function<double(double, double)> func);
    /**
     * @brief Метод построения массива ребер
     * @return Массив ребер
     */
    std::vector<Segment> evalEdges();
    /**
     * @brief Метод вычисления оптимальных координат узла при помощи минимизации функционала експоненты площади с учтом знака
     * @param i Номер (код) узла
     * @param t Параметр экспоненты площади (по умолчанию -10)
     * @return Координаты оптимального положения узлы
     */
    Point2D evalOptimalPosition(const UInteger &i, double t=3.0);
    void optimizeBorder(std::function<double(double, double)> func, int iiter = 4, double level = 0.0);
protected:
    std::vector<Triangle> element_; //!< Массив элементов
    typedef std::vector<Triangle>::iterator ElementIterator;
};
}

#endif // TRIANGLEMESH2D_H
