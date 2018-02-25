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
#include "trianglemesh2d.h"
#include "tetrahedralmesh3d.h"

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
     * @brief Конструктор создает копию объекта, переданого по указателю
     * @param mesh Указатель на объект для копирования
     */
    TriangleMesh3D(const TriangleMesh3D *mesh);
    /**
     * @brief Метод для создания равномерной сетки в цилиндрических координатах
     * @param rCount Количество элементов вдоль радиуса
     * @param lCount Количество элементов вдоль образующей
     * @param radius Радиус
     * @param length Высота
     */
    void cylinderDomain(const UInteger &rCount, const UInteger &lCount, const double &radius, const double &length);
    /**
     * @brief Метод для создания равномерной сетки в цилиндрических координатах с адаптацией под область, ограниченную функцией
     * @param rCount Количество элементов вдоль радиуса
     * @param lCount Количество элементов вдоль образующей
     * @param radius Радиус
     * @param length Высота
     * @param func Функция, положительная во внутренних точка и отрицательная во внешних
     */
    void cylinderDomain(const UInteger &rCount, const UInteger &lCount, const double &radius, const double &length, std::function<double(double, double, double)> func);
    /**
     * @brief Метод для создания равномерной сетки в конических координатах
     * @param rCount Количество элементов вдоль радиуса
     * @param lCount Количество элементов вдоль образующей
     * @param bottom_radius Нижний радиус
     * @param top_radius Верхний радиус
     * @param length Высота
     */
    void coneDomain(const UInteger &rCount, const UInteger &lCount, const double &bottom_radius, const double &top_radius, const double &length);
    /**
     * @brief Метод для создания равномерной сетки в конических координатах с адаптацией под область, ограниченную функцией
     * @param rCount Количество элементов вдоль радиуса
     * @param lCount Количество элементов вдоль образующей
     * @param bottom_radius Нижний радиус
     * @param top_radius Верхний радиус
     * @param length Высота
     * @param func Функция, положительная во внутренних точка и отрицательная во внешних
     */
    void coneDomain(const UInteger &rCount, const UInteger &lCount, const double &bottom_radius, const double &top_radius, const double &length, std::function<double(double, double, double)> func);
    void parametricDomain(const UInteger &uCount, const UInteger &vCount, std::function<Point3D(double, double)> domainFunction, std::function<double(double, double, double)> rfunc);
    /**
     * @brief Метод марширующих кубиков
     * @param xCount Количество точек дискретизации по оси абсцисс
     * @param yCount Количество точек дискретизации по оси ординат
     * @param zCount Количество точек дискретизации по оси аппликат
     * @param xMin Абсцисса минимального угла
     * @param yMin Ордината минимального угла
     * @param zMin Аппликата минимального угла
     * @param width Ширина
     * @param height Высота
     * @param depth Глубина
     * @param func Функция, положительная во внутренних точка и отрицательная во внешних
     * @param level Уровень поверхности
     */
    void marchingCubes(const UInteger &xCount, const UInteger &yCount, const UInteger &zCount, const double &xMin, const double &yMin, const double &zMin, const double &width, const double &height, const double &depth, std::function<double (double, double, double)> func, double level = 0.0, bool slice_x = false, bool slice_y = false, bool slice_z = false);
    /**
     * @brief Метод марширующих тетраэдров
     * @param xCount Количество точек дискретизации по оси абсцисс
     * @param yCount Количество точек дискретизации по оси ординат
     * @param zCount Количество точек дискретизации по оси аппликат
     * @param xMin Абсцисса минимального угла
     * @param yMin Ордината минимального угла
     * @param zMin Аппликата минимального угла
     * @param width Ширина
     * @param height Высота
     * @param depth Глубина
     * @param func Функция, положительная во внутренних точка и отрицательная во внешних
     * @param level Уровень поверхности
     */
    void marchingTetrahedrons(const UInteger &xCount, const UInteger &yCount, const UInteger &zCount, const double &xMin, const double &yMin, const double &zMin, const double &width, const double &height, const double &depth, std::function<double (double, double, double)> func, double level = 0.0, bool slice_x = false, bool slice_y = false, bool slice_z = false);
    /**
     * @brief Метод построения стеки с использованием фоновых тетраэдров
     * @param mesh Сетка тетраэдров
     * @param func Функция области
     * @param level Уровень нуля
     * @param smooth Количество итераций сглаживания
     * @param optimize Количество итераций оптимизации
     */
    void backgroundGrid(const TetrahedralMesh3D *mesh, std::function<double(double, double, double)> func, double level = 0.0, int smooth = 0, int optimize = 0, bool useFlip = true);
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
     * @brief Добавить элемент к сетке
     * @param Массив ссылок (номеров) на узлы
     */
    void addElement(const std::vector<UInteger> &nodes_ref);
    void addElementOrdered(const Point3D &t0, const Point3D &t1, const Point3D &t2, double epsilon = epsilon_);
    /**
     * @brief area Вычислить площадь элемента
     * @param number Номер элемента
     * @return Площадь элемента
     */
    virtual double area(const UInteger &number) const;
    /**
     * @brief Вычислить площадь треугольника, заданного тремя вершинами (формула Герона)
     * @param p0 Первая вершина
     * @param p1 Вторая вершина
     * @param p2 Третья вершина
     * @return Площадь треугольника
     */
    virtual double area(const Point3D &p0, const Point3D &p1, const Point3D &p2) const;
    /**
     * @brief Вычислить значение минимального угла в элементе
     * @param elNum Номер элемента
     * @return Минимальный угол элемента (радианы)
     */
    double minAngle(const UInteger &elNum);
    /**
     * @brief Метод очищает информацию об елементах
     */
    virtual void clearElements();
    /**
     * @brief Операция добавления сетки к существующей на основе элементарного объединения
     * @param mesh Указатель на сетку для добавления
     */
    void add(const TriangleMesh3D *mesh);
    /**
     * @brief Функция принадлежности точки контуру
     * @param x Абсцисса точки
     * @param y Ордината точки
     * @param z Аппликата точки
     * @return Расстояние до ближайшей точки границы со знаком "+" для внутренних точек
     */
    double cfunction(const double &x, const double &y, const double &z);
    /**
     * @brief Процедура сглаживания Лапласа
     * @param func Функция области
     * @param level Линия уровня
     * @param iter_num Количесво итераций
     */
    void laplacianSmoothing(std::function<double(double, double, double)> func, double level = 0, int iter_num = 1, bool useFlip = true);
    /**
     * @brief Процедура сглаживания на основе минимизации функционала расстояния-длины
     * @param func Функция области
     * @param level Линия уровня
     * @param iter_num Количесво итераций
     */
    void distlenSmoothing(std::function<double(double, double, double)> func, double level = 0, int iter_num = 1, bool useFlip = true);

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
    /**
     * @brief Метод проверяет попадание точки P в сферу, для которой описанная окружность треугольника ABC является диаметральным сечением
     * @param P Координаты точки для проверки
     * @param A Первая точка треугольника
     * @param B Втораяточка треуольника
     * @param C Третья точка треугольника
     * @return true, если точка попадает в сферу
     */
    bool inCircumSphere(const Point3D &P, const Point3D &A, const Point3D &B, const Point3D &C);
    bool inCircumCylinder(const Point2D &P, const Point2D &A, const Point2D &B, const Point2D &C, std::function<Point3D(double, double)> domainFunction);
    Point2D circumCylinderCenter(const Point2D &A, const Point2D &B, const Point2D &C, std::function<Point3D(double, double)> domainFunction);
    bool insertDelaunayNode(const Point2D &point, const NodeType &type, TriangleMesh2D::Triangulation &triangulation, std::function<Point3D(double, double)> domainFunction);
    void flip();
    typedef struct
    {
        Point3D a, b, c;
    } CooTriangle;
    /**
     * @brief Процедура удаления из списка треугольников, которые со стороной меньше delta
     * @param cootriangles Список треугольников, заданных своими координатами
     * @param func Функция области
     * @param level Уровень границы
     * @param delta Минимально допустимый размер стороны
     */
    void clearCooTriangles(std::list<CooTriangle> &cootriangles, std::function<double (double, double, double)> func, double level = 0.0, double delta = epsilon_);
protected:
    std::vector<Triangle> element_; //!< Массив элементов
};
}
#endif // TRIANGLEMESH3D_H
