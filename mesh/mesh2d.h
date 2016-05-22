/**
  * @author Сергей Чопоров
  * @date 11/02/2014
  * @version 1.0.2
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef MESH2D_H
#define MESH2D_H

#include <vector>
#include <functional>

#include "mesh.h"
#include "node2d.h"

namespace msh
{
/**
 * @brief Класс Mesh2D - неструтурированная сетка элементов.
 * Неструктурированная сетка представляется двумя списками: вершин и связей (элементов)
 * Вспомогательным является список множест смежных в узле элементов
 * В данном классе не определена структура элементов, следовательно, отсутсвует их список
 * @see Mesh, Node2D, Point2D
 */
class Mesh2D : public Mesh
{
public:
    Mesh2D(const Mesh2D *mesh);
    /**
     * @brief Количество узлов
     * @return Количество узлов в сетке
     */
    virtual UInteger nodesCount() const;
    /**
     * @brief Координаты узла сетки
     * @param number Номер узла сетки
     * @return Указатель на координаты узла сетки с номером number
     */
    virtual PointPointer node(const UInteger &number) const;
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
     * @return Размерность пространства, в котором определена сетка
     */
    virtual int dimesion() const;
    /**
     * @brief Тип узла сетки
     * @param number Номер узла сетки
     * @return Тип узла сетки с номером number
     */
    virtual NodeType nodeType(const UInteger &number) const;
    /**
     * @brief Метод очищает информацию об узлах сетки
     */
    virtual void clearNodes();
    /**
     * @brief Отразить сетку по вертикали
     */
    void flipVertically();
    /**
     * @brief Отразить сетку по горизонтали
     */
    void flipHorizontally();
    /**
     * @brief Зеркально отразить сетку по вертикали
     */
    void mirrorVertically();
    /**
     * @brief Зеркально отразить сетку по горизонтали
     */
    void mirrorHorizontally();
    /**
     * @brief Изменить направление обхода во всех элементах сетки
     */
    virtual void directionChange() = 0;
    /**
     * @brief Подсчитать площадь дискретной модели
     * @return Сумма площадей конечных элементов
     */
    virtual double area() const;
    /**
     * @brief area Подсчитать площадь конечного элемента
     * @param number Номер элемента
     * @return Площадь элемента
     */
    virtual double area(const UInteger &number) const = 0;
    /**
     * @brief Добавить узел, заданный точкой, без проверки на наличие в массиве узлов
     * @param point Координаты узла для вставки
     * @param type Тип узла
     * @return Номер узла в массиве узлов
     */
    UInteger pushNode(const Point2D &point, const NodeType &type);
    /**
     * @brief Добавить узел в сетку
     * @param point Указатель на координаты
     * @param type Тип узла
     * @return Номер вставленного узла в сетке
     */
    UInteger pushNode(PointPointer point, const NodeType &type);
    /**
     * @brief Добавить узел, заданный точкой, с использованием проверки на наличие узла с такими координатами
     * @param point Координаты узла для вставки
     * @param type Тип узла
     * @return Номер узла в массиве узлов
     */
    UInteger addNode(const Point2D &point, const NodeType &type, double epsilon = epsilon_, std::function<double(Point2D, Point2D)> distance = nullptr);
    /**
     * @brief Обновить параметры области определения сетки (xMin, xMax, yMin, yMax)
     */
    virtual void updateDomain();
    /**
     * @brief Количество элементов, смежных в узле с заданным номером
     * @param nodeNumber Номер узла
     * @return Количество элементов, смежных в узле
     */
    UInteger adjacentCount(const UInteger &nodeNumber) const;
    /**
     * @brief Метод возвращает множество номеров смежных в узле элементов
     * @param nodeNumber Номер узла
     * @return Множество смежных элементов (их номера)
     */
    AdjacentSet adjacent(const UInteger &nodeNumber) const;
    /**
     * @brief Метод возвращает двумерные координаты узла
     * @param number Номер узла
     * @return Двумерные координаты узла
     */
    Point2D point2d(const UInteger &number) const;
    /**
     * @brief Метод изменяет координаты узла на заданные
     * @param number Номер узла
     * @param p Координаты
     */
    void setPoint(const UInteger &number, const Point2D &p);
    /**
     * @brief Двоичный поиск границы функции, пересекаемой отрезком
     * @param p0 Первая точка отрезка, пересекающего границу
     * @param p1 Вторая точка отрезка, пересекающего границу
     * @param func Функция области
     * @return Координаты пересечения отрезка и границы области
     */
    static Point2D binary(Point2D p0, Point2D p1, std::function<double(double, double)> func);
protected:
    std::vector<Node2D> node_; //!< массив узлов
    double xMin_; //!< минимальное значение ординаты
    double xMax_; //!< максимальное значение ординаты
    double yMin_; //!< минимальное значение абсциссы
    double yMax_; //!< максимальное значение абсциссы
};
}

#endif // MESH2D_H
