/**
 * @author Сергей Чопоров
 * @version 1.0.4
 * @date 19.11.2012
 * @copyright Copyright 2012 Sergey Choporov. All rights reserved.
 * This project is released under the GNU Lesser General Public License.
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
     * @param node_0 Номер (код) нулевого узла
     * @param node_1 Номер (код) первого узла
     * @param node_2 Номер (код) второго узла
     * @param node_3 Номер (код) третьего узла
     */
    Quadrilateral(const UInteger &node_0, const UInteger &node_1, const UInteger &node_2, const UInteger &node_3);
    /**
     * @brief Конструктор копирования
     * @param quad Объект-четырехугольник для копирования
     */
    Quadrilateral(const Quadrilateral &quad);
    /**
     * @brief Оператор для доступа к номеру (коду) вершины по ее индексу в четырхеугольнике
     * @param i Индекс вершины четырехугольника
     * @return Номер (код) узла, на который ссылается вершина
     */
    UInteger &operator [] (int i);
    /**
     * @brief Оператор (коснтантный) для доступа к номеру (коду) вершины по ее индексу в четырхеугольнике
     * @param i Индекс вершины четырехугольника
     * @return Номер (код) узла, на который ссылается вершина
     */
    const UInteger &operator [] (int i) const;
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
     * @brief Номер узла (код), на который ссылается вершина элемента
     * @param i Номер вершины
     * @return Номер узла (код), на который ссылается вершина
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
    UInteger vertex_[4]; //!< Вершины четырехугольника (номера, хеш-коды)
};
}

#endif // QUADRILATERAL_H
