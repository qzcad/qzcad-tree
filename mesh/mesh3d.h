/**
  * @author Сергей Чопоров
  * @date 19/06/2014
  * @version 1.0.1
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef MESH3D_H
#define MESH3D_H

#include <vector>
#include <functional>

#include "mesh.h"
#include "node3d.h"

namespace msh {
/**
 * @brief Класс Mesh3D - неструтурированная сетка элементов.
 * Неструктурированная сетка представляется двумя списками: вершин и связей (элементов)
 * Вспомогательным является список множест смежных в узле элементов
 * В данном классе не определена структура элементов, следовательно, отсутсвует их список
 * @see Mesh, Node3D, Point3D
 */
class Mesh3D: public Mesh
{
public:
    Mesh3D(const Mesh3D *mesh);
    /**
     * @brief Количество узлов
     * @return Количество узлов в сетке
     */
    virtual UInteger nodesCount() const;
    /**
     * @brief Координаты узла сетки
     * @param number Номер узла сетки
     * @return Указатель на копию координат узла сетки с номером number
     */
    virtual PointPointer node(const UInteger &number) const;
    /**
     * @brief Тип узла сетки
     * @param number Номер узла сетки
     * @return Тип узла сетки с номером number
     */
    virtual NodeType nodeType(const UInteger &number) const;
    /**
     * @brief Минимальное значение ординаты точек сетки
     * @return Минимальное значение ординаты точек сетки
     */
    virtual double xMin() const;
    /**
     * @brief Максимальное значение ординаты точек сетки
     * @return Максимальное значение ординаты точек сетки
     */
    virtual double xMax() const;
    /**
     * @brief Минимальное значение абсциссы точек сетки
     * @return Минимальное значение абсциссы точек сетки
     */
    virtual double yMin() const;
    /**
     * @brief Максимальное значение абсциссы точек сетки
     * @return Максимальное значение абсциссы точек сетки
     */
    virtual double yMax() const;
    /**
     * @brief Минимальное значение аппликаты точек сетки
     * @return Минимальное значение аппликаты точек сетки
     */
    virtual double zMin() const;
    /**
     * @brief Максимальное значение аппликаты точек сетки
     * @return Максимальное значение аппликаты точек сетки
     */
    virtual double zMax() const;
    /**
     * @brief Размерность пространства
     * @return Размерность пространства, в котором определена сетка (3)
     */
    virtual int dimesion() const;
    /**
     * @brief Метод очищает информацию об узлах сетки
     */
    virtual void clearNodes();
    /**
     * @brief Вычислить площадь поверхности дискретной модели
     * @return Площадь поверхности дискретной модели
     */
    virtual double surfaceArea() const = 0;
    /**
     * @brief Добавить узел, заданный точкой, без проверки на наличие в массиве узлов
     * @param point Координаты узла для вставки
     * @param type Тип узла
     * @return Номер узла в массиве узлов
     */
    UInteger pushNode(const Point3D &point, const NodeType &type);
    /**
     * @brief Добавить узел в сетку
     * @param point Указатель на координаты
     * @param type Тип узла
     * @return Номер вставленного узла в сетке
     */
    UInteger pushNode(PointPointer point, const NodeType &type);
    /**
     * @brief Добавить узел, заданный точкой, с проверкой на наличие в массиве узлов
     * @param point Координаты узла для вставки
     * @param type Тип узла
     * @return Номер узла в массиве узлов
     */
    UInteger addNode(const Point3D &point, const NodeType &type, double epsilon = epsilon_);
    /**
     * @brief Добавить узел с проверкой на наличие в массиве узлов
     * @param node Коордианты и тип узла
     * @return Номер узла в массиве узлов
     */
    UInteger addNode(const Node3D &node, double epsilon = epsilon_);
    /**
     * @brief Обновить параметры области определения сетки (xMin, xMax, yMin, yMax, zMin, zMax)
     */
    virtual void updateDomain();
    /**
     * @brief Метод возвращает множество номеров смежных в узле элементов
     * @param nodeNumber Номер узла
     * @return Множество смежных элементов (их номера)
     */
    virtual AdjacentSet adjacent(const UInteger &nodeNumber) const;
    /**
     * @brief Количество элементов, смежных в узле с заданным номером
     * @param nodeNumber Номер узла
     * @return Количество элементов, смежных в узле
     */
    virtual UInteger adjacentCount(const UInteger &nodeNumber) const;
    /**
     * @brief Двоичный поиск границы функции, пересекаемой отрезком
     * @param p0 Первая точка отрезка, пересекающего границу
     * @param p1 Вторая точка отрезка, пересекающего границу
     * @param func Функция области
     * @param level Линия уровня
     * @return Координаты пересечения отрезка и границы области
     */
    static Point3D binary(Point3D p0, Point3D p1, std::function<double(double, double, double)> func, double level = 0.0);
    /**
     * @brief Поиск ближайшей к точке point граничной при помощи градиента, вычисленного разностной схемой с шагом h
     * @param point Точка, для которой необходимо найти ближайшую граничную
     * @param func Функция области
     * @param h Шаг для вычисления частных производных
     * @param level Линия уровня границы
     * @return Координаты граничной точки
     */
    static Point3D findBorder(Point3D point, std::function<double (double, double, double)> func, double h, double level = 0.0);
    /**
     * @brief Поиск ближайшей граничной точки к точке, определенной l-координатами треугольника: point = ((1.0 - lb - lc) * a) + (lb * b) + (lc * c)
     * @param a Первая вершина треугольника (в порядке обхода)
     * @param b Вторая вершина треугольника (в порядке обхода)
     * @param c Третья вершина треугольника (в порядке обхода)
     * @param func Функция области
     * @param lb Первая l-координата
     * @param lc Вторая l-координата
     * @param level Линия уровня
     * @return Координаты точки на границе области
     */
    static Point3D findBorder(const Point3D &a, const Point3D &b, const Point3D &c, std::function<double (double, double, double)> func, double lb = 0.33333, double lc = 0.33333, double level = 0.0);
    /**
     * @brief Вычисление расстояния до границы от точки, определенной l-координатами треугольника: point = ((1.0 - lb - lc) * a) + (lb * b) + (lc * c)
     * @param a Первая вершина треугольника (в порядке обхода)
     * @param b Вторая вершина треугольника (в порядке обхода)
     * @param c Третья вершина треугольника (в порядке обхода)
     * @param func Функция области
     * @param lb Первая l-координата
     * @param lc Вторая l-координата
     * @param level Линия уровня
     * @return Расстояние до границы области
     */
    static double distToBorder(const Point3D &a, const Point3D &b, const Point3D &c, std::function<double (double, double, double)> func, double lb = 0.33333, double lc = 0.33333, double level = 0.0);
    /**
     * @brief Вычислить градиент в точке (четырехточечная схема)
     * @param func Функция двух переменных
     * @param p Точка, в которой вычисляется градиент
     * @param h Шаг разностной схемы
     * @return Координаты градиента
     */
    static Point3D grad(std::function<double(double, double, double)> func, const Point3D &p, const double &h);
    /**
     * @brief Метод для построения вектора значений заданной функции в узлах сетки
     * @param func Указатель функции трех переменных
     */
    void evalNodalValues(std::function<double (double, double, double)> func);
protected:
    std::vector<Node3D> node_; //!< массив узлов
    double xMin_; //!< минимальное значение ординаты
    double xMax_; //!< максимальное значение ординаты
    double yMin_; //!< минимальное значение абсциссы
    double yMax_; //!< максимальное значение абсциссы
    double zMin_; //!< минимальное значение аппликаты
    double zMax_; //!< максимальное значение аппликаты
};
}

#endif // MESH3D_H
