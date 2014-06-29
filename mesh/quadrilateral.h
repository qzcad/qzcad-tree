/**
 * @author Сергей Чопоров
 * @version 1.0.4
 * @date 19.11.2012
 */
#ifndef QUADRILATERAL_H
#define QUADRILATERAL_H

#include "element.h"

namespace msh
{
/**
 * @brief Класс Quadrilateral - четырехугольник
 * Обход узлов против часовой стрелки
 * vertex[3] vertex[2]
 *  *-------*
 *  |       |
 *  |       |
 *  *-------*
 * vertex[0] vertex[1]
 */
class Quadrilateral: public Element
{
public:
    /**
     * @brief Конструктор
     * @param node_0 Номер нулевого узла
     * @param node_1 Номер первого узла
     * @param node_2 Номер второго узла
     * @param node_3 Номер третьего узла
     */
    Quadrilateral(const UInteger &node_0, const UInteger &node_1, const UInteger &node_2, const UInteger &node_3);
    /**
     * @brief Конструктор копирования
     * @param quad Объект-четырехугольник для копирования
     */
    Quadrilateral(const Quadrilateral &quad);
    /**
     * @brief operator []
     * @param i Индекс вершины четырехугольника
     * @return Номер узла, на который ссылается вершина
     */
    UInteger &operator [] (int i);
    /**
     * @brief operator =
     * @param quad Объект-четырехугольник для копирования
     * @return Копию четырехугольника
     */
    Quadrilateral &operator = (const Quadrilateral &quad);
    /**
     * @brief nodesCount Количество вершин в элементе
     * @return 4
     */
    virtual int verticesCount() const;
    /**
     * @brief Номер узла, на который ссылается вершина элемента
     * @param i Номер вершины
     * @return Номер узла, на который ссылается вершина
     */
    virtual UInteger vertexNode(int i) const;
    /**
     * @brief Количество граней, соответвующих элементу (необходимо для визуализации)
     * @return Количество граней (1)
     */
    virtual int facesCount() const;
    /**
     * @brief face Грань элемента
     * @param i Номер грани
     * @return Грань элемента, номер которойс равен i
     */
    virtual UIntegerVector face(const int &i) const;
private:
    UInteger vertex_[4]; //!< Вершины четырехугольника
};
}

#endif // QUADRILATERAL_H
