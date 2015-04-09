/**
 * @author Сергей Чопоров
 * @version 1.0.0
 * @date 9.04.2015
 * @copyright Copyright 2015 Sergey Choporov. All rights reserved.
 * This project is released under the GNU Lesser General Public License.
 */
#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "element.h"

namespace msh
{
/**
 * @brief Класс Triangle - треугольник.
 * Обход узлов против часовой стрелки.
 *      vertex[2]
 *         *
 *        / \
 *       /   \
 *      /     \
 *     *-------*
 * vertex[0]  vertex[1]
 */
class Triangle : public Element
{
public:
    /**
     * @brief Конструктор
     * @param node_0 Номер (код) нулевого узла
     * @param node_1 Номер (код) первого узла
     * @param node_2 Номер (код) второго узла
     */
    Triangle(const UInteger &node_0, const UInteger &node_1, const UInteger &node_2);
    /**
     * @brief Конструктор копирования
     * @param triangle Экземпляр объекта для копирования
     */
    Triangle(const Triangle &triangle);
    /**
     * @brief Оператор для доступа к номеру (коду) вершины по ее индексу в треугольнике
     * @param i Индекс вершины треугольника
     * @return Номер (код) узла, на который ссылается вершина
     */
    UInteger &operator [] (int i);
    /**
     * @brief Оператор (коснтантный) для доступа к номеру (коду) вершины по ее индексу в треугольнике
     * @param i Индекс вершины треугольника
     * @return Номер (код) узла, на который ссылается вершина
     */
    const UInteger &operator [] (int i) const;
    /**
     * @brief Оператор присваивания
     * @param triangle Объект-треугольник для копирования
     * @return Копию треугольника
     */
    Triangle &operator = (const Triangle &triangle);
    /**
     * @brief nodesCount Количество вершин в элементе
     * @return 3
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
    UInteger vertex_[3]; //!< Вершины треугольника (номера или хеш-коды)
};
}

#endif // TRIANGLE_H
