  /**
  * @author Сергей Чопоров
  * @date 15/10/2012
  * @version 1.0.6
  */
#ifndef MESH_H
#define MESH_H

#include "pointpointer.h"
#include "elementpointer.h"
#include "nodetype.h"
namespace msh
{
/**
 * @brief Класс Mesh - абстракция сетки элементов
 * В классе Mesh не определены тип сетки и форма элемента
 * @see Element, Point, PointPointer, ElementPointer, UInteger, Floating, NodeType
 */
class Mesh
{

public:
    /**
     * @brief Количество узлов
     * @return Количество узлов в сетке
     */
    virtual UInteger nodesCount() const = 0;
    /**
     * @brief Координаты узла сетки
     * @param number Номер узла сетки
     * @return Указатель на копию координат узла сетки с номером number
     */
    virtual PointPointer node(const UInteger &number) const = 0;
    /**
     * @brief Тип узла сетки
     * @param number Номер узла сетки
     * @return Тип узла сетки с номером number
     */
    virtual NodeType nodeType(const UInteger &number) const = 0;
    /**
     * @brief Количество элементов
     * @return Количество элементов в сетке
     */
    virtual UInteger elementsCount() const = 0;
    /**
     * @brief Элемет сетки
     * @param number Номер элемента сетки
     * @return Указатель на копию элемента сетки
     */
    virtual ElementPointer element(const UInteger &number) const = 0;
    /**
     * @brief Минимальное значение ординаты точек сетки
     * @return Минимальное значение ординаты точек сетки
     */
    virtual Floating xMin() const = 0;
    /**
     * @brief Максимальное значение ординаты точек сетки
     * @return Максимальное значение ординаты точек сетки
     */
    virtual Floating xMax() const = 0;
    /**
     * @brief Минимальное значение абсциссы точек сетки
     * @return Минимальное значение абсциссы точек сетки
     */
    virtual Floating yMin() const = 0;
    /**
     * @brief Максимальное значение абсциссы точек сетки
     * @return Максимальное значение абсциссы точек сетки
     */
    virtual Floating yMax() const = 0;
    /**
     * @brief Минимальное значение аппликаты точек сетки
     * @return Минимальное значение аппликаты точек сетки
     */
    virtual Floating zMin() const = 0;
    /**
     * @brief Максимальное значение аппликаты точек сетки
     * @return Максимальное значение аппликаты точек сетки
     */
    virtual Floating zMax() const = 0;
    /**
     * @brief Размерность пространства
     * @return Размерность пространства, в котором определена сетка
     */
    virtual int dimesion() const = 0;
    /**
     * @brief Значение некоторой функции, определенной на элементе
     * @param number Номер элемента
     * @return Значение, соответствующее элементу
     */
    virtual Floating elementValue(const UInteger &number) const = 0;
    /**
     * @brief Значение некоторой функции, определенной на зле
     * @param number Номер узла
     * @return Значение, соответствующее узлу
     */
    virtual Floating nodeValue(const UInteger &number) const = 0;

    virtual bool isBorderElement(const UInteger &number) const = 0;
};
}
#endif // MESH_H
