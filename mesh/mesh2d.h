/**
  * @author Сергей Чопоров
  * @date 11/02/2014
  * @version 1.0.2
  */
#ifndef MESH2D_H
#define MESH2D_H

#include <vector>

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
    Mesh2D();
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
    virtual Floating xMin() const;
    /**
     * @brief Максимальное значение ординаты точек сетки
     * @return Максимальное значение ординаты точек сетки
     */
    virtual Floating xMax() const;
    /**
     * @brief Минимальное значение абсциссы точек сетки
     * @return Минимальное значение абсциссы точек сетки
     */
    virtual Floating yMin() const;
    /**
     * @brief Максимальное значение абсциссы точек сетки
     * @return Максимальное значение абсциссы точек сетки
     */
    virtual Floating yMax() const;
    /**
     * @brief Минимальное значение аппликаты точек сетки
     * @return Минимальное значение аппликаты точек сетки
     */
    virtual Floating zMin() const;
    /**
     * @brief Максимальное значение аппликаты точек сетки
     * @return Максимальное значение аппликаты точек сетки
     */
    virtual Floating zMax() const;
    /**
     * @brief Размерность пространства
     * @return Размерность пространства, в котором определена сетка
     */
    virtual int dimesion() const;
    /**
     * @brief Значение некоторой функции, опеределенной на элемента
     * @param number Номер элемента
     * @return Значение равное number
     */
    virtual Floating elementValue(const UInteger &number) const;
    /**
     * @brief Значение некоторой функции, определенной на зле
     * @param number Номер узла
     * @return Значение, соответствующее узлу
     */
    virtual Floating nodeValue(const UInteger &number) const;
    /**
     * @brief Тип узла сетки
     * @param number Номер узла сетки
     * @return Тип узла сетки с номером number
     */
    virtual NodeType nodeType(const UInteger &number) const;
protected:
    /**
     * @brief Добавить узел, заданный точкой, без проверки на наличие в массиве узлов
     * @param point Координаты узла для вставки
     * @param type Тип узла
     * @return Номер узла в массиве узлов
     */
    UInteger pushNode(const Point2D &point, const NodeType &type);
    /**
     * @brief Добавить узел, заданный точкой, с использованием проверки на наличие узла с такими координатами
     * @param point Координаты узла для вставки
     * @param type Тип узла
     * @return Номер узла в массиве узлов
     */
    UInteger addNode(const Point2D &point, const NodeType &type);
protected:
    std::vector<Node2D> node_; //!< массив узлов
    Floating xMin_; //!< минимальное значение ординаты
    Floating xMax_; //!< максимальное значение ординаты
    Floating yMin_; //!< минимальное значение абсциссы
    Floating yMax_; //!< максимальное значение абсциссы
};
}

#endif // MESH2D_H
