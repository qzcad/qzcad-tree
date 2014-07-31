/**
  * @author Сергей Чопоров
  * @date 19/06/2014
  * @version 1.0.1
  */
#ifndef MESH3D_H
#define MESH3D_H

#include <vector>

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
    Mesh3D();
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
     * @return Размерность пространства, в котором определена сетка (3)
     */
    virtual int dimesion() const;
    /**
     * @brief Значение некоторой функции, определенной на элементе
     * @param number Номер элемента
     * @return Значение, соответствующее элементу
     */
    virtual Floating elementValue(const UInteger &number) const;
    /**
     * @brief Значение некоторой функции, определенной на зле
     * @param number Номер узла
     * @return Значение, соответствующее узлу
     */
    virtual Floating nodeValue(const UInteger &number) const;
    /**
     * @brief Вычислить площадь поверхности дискретной модели
     * @return Площадь поверхности дискретной модели
     */
    virtual Floating surfaceArea() const = 0;
    /**
     * @brief Добавить узел, заданный точкой, без проверки на наличие в массиве узлов
     * @param point Координаты узла для вставки
     * @param type Тип узла
     * @return Номер узла в массиве узлов
     */
    UInteger pushNode(const Point3D &point, const NodeType &type);
    /**
     * @brief Добавить узел, заданный точкой, с проверкой на наличие в массиве узлов
     * @param point Координаты узла для вставки
     * @param type Тип узла
     * @return Номер узла в массиве узлов
     */
    UInteger addNode(const Point3D &point, const NodeType &type);
protected:
    std::vector<Node3D> node_; //!< массив узлов
    Floating xMin_; //!< минимальное значение ординаты
    Floating xMax_; //!< максимальное значение ординаты
    Floating yMin_; //!< минимальное значение абсциссы
    Floating yMax_; //!< максимальное значение абсциссы
    Floating zMin_; //!< минимальное значение аппликаты
    Floating zMax_; //!< максимальное значение аппликаты
};
}

#endif // MESH3D_H
