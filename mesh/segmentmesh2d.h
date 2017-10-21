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
     * @brief Метод создает равномерную дискретную модель заданной линии уровня на основе сегментов для области, определенной функционально
     * @param xCount Количество узлов вдоль оси абсцисс
     * @param yCount Количество узлов вдоль оси ординат
     * @param xMin Абсцисса нижнего левого угла прямоугольной области
     * @param yMin Ордината нижнего левого угла прямоугольной области
     * @param width Ширина прямоугольной области
     * @param height Высота прямоугольной области
     * @param func Функция области
     * @param charPoint Список характерных точек
     * @param level Линия уровня для дискретизации
     * @param isOptimized Если true, то сетка сгущается на участках с наибольшей кривизной
     * @param distance Указатель на фукнцию растояния (если дискретизация не в двумерной декартовой системе)
     */
    void functionalDomain(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double(double, double)> func, std::list<Point2D> charPoint, double level = 0.0, bool isOptimized = true, std::function<double(Point2D, Point2D)> distance = nullptr);
    /**
     * @brief Метод создает равномерную дискретную модель на основе сегментов для области, определенной функционально (случай контакта двух тел)
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
     */
    void functionalDomain(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double(double, double)> func_a, std::function<double(double, double)> func_b, std::list<Point2D> charPoint, double delta = -1.0);
    /**
     * @brief Метод создает дискретную модель заданного числа линий уровня R-функция
     * @param xCount Количество узлов вдоль оси абсцисс
     * @param yCount Количество узлов вдоль оси ординат
     * @param xMin Абсцисса нижнего левого угла прямоугольной области
     * @param yMin Ордината нижнего левого угла прямоугольной области
     * @param width Ширина прямоугольной области
     * @param height Высота прямоугольной области
     * @param func Функция области
     * @param charPoint Список характерных точек (уровня 0.0)
     * @param contours Количество линий уровня
     * @param isOptimized Если true, то сетка сгущается на участках с наибольшей кривизной
     */
    void contourGraph(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double(double, double)> func, std::list<Point2D> charPoint, int contours = 4, bool isOptimized = true);
    /**
     * @brief Метод создает дискретную модель заданного числа внутренних фронтов границы R-функция
     * @param xCount Количество узлов вдоль оси абсцисс
     * @param yCount Количество узлов вдоль оси ординат
     * @param xMin Абсцисса нижнего левого угла прямоугольной области
     * @param yMin Ордината нижнего левого угла прямоугольной области
     * @param width Ширина прямоугольной области
     * @param height Высота прямоугольной области
     * @param func Функция области
     * @param charPoint Список характерных точек (уровня 0.0)
     * @param contours Количество линий уровня
     * @param isOptimized Если true, то сетка сгущается на участках с наибольшей кривизной
     */
    void frontGraph(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, std::function<double(double, double)> func, std::list<Point2D> charPoint, int contours = 4, bool isOptimized = true);
    void parametricDomain(const UInteger &count, const double &tmin, const double &tmax, std::function<Point2D (double)> domainFunction);
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
     * @param segment Экземляр элемента для добавления в сетку
     */
    void addElement(const Segment &segment);
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
    /**
     * @brief Метод очищает информацию об елементах
     */
    virtual void clearElements();
    /**
     * @brief Метод проверяет пересечение заданного отрезка сегментом сетки (! строго его внутренней частью)
     * @param p0 Первая точка отрезка
     * @param p1 Вторая точка отрезка
     * @param number Номер первого сегмента, с которым найдено пересечение
     * @return true, если был найден походящий сегмент
     */
    bool isCrossedElement(const Point2D &p0, const Point2D &p1, UInteger &number);
    /**
     * @brief Метод репроекции средины сегмента на границу
     * @param number Номер сегмента (элемента)
     * @param func Функция области
     * @return Координаты нового узла
     */
    Point2D refineMidpoint(const UInteger &number, std::function<double(double, double)> func);
    /**
     * @brief Метод проверяет попадание точки в окружность, радиусом которой приходится один из сегментов (элементов)
     * @param point Координаты точки
     * @param number Номер первого элемента, для которого обнаружено вхождение
     * @return true, если был найден походящий сегмент
     */
    bool isEncroached(const Point2D &point, UInteger &number);
    /**
     * @brief Процедура сглаживания Лапласа (применяется только к узла нулевой линии уровня)
     * @param iter_num Количесво итераций
     */
    void laplacianSmoothing(std::function<double(double, double)> func, double level = 0.0, int iter_num = 4);
    void distlenSmoothing(std::function<double(double, double)> func, double level = 0.0, int iter_num = 8);
    /**
     * @brief Функция принадлежности точки контуру
     * @param x Абсцисса точки
     * @param y Ордината точки
     * @return Расстояние до ближайшей точки границы со знаком "+" для внутренних точек
     */
    double cfunction(const double &x, const double &y);
protected:
    void cellContours(const Point2D &p0, const Point2D &p1, const Point2D &p2, const Point2D &p3, const double &v0, const double &v1, const double &v2, const double &v3, std::function<double(double, double)> func, double level = 0.0, std::function<double(Point2D, Point2D)> distance = nullptr);
private:
    std::vector<Segment> element_; //!< Массив элементов
    typedef std::vector<Segment>::iterator ElementIterator;
};
}

#endif // SEGMENTMESH2D_H
