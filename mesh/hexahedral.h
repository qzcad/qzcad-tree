/**
  * @author Сергей Чопоров
  * @date 19/06/2014
  * @version 1.0.1
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef HEXAHEDRAL_H
#define HEXAHEDRAL_H

#include "element.h"

namespace msh {
/**
 * @brief Класс Hexahedral - шестигранный элемент в пространстве
 * ^
 * | Y
 *       4___________7
 *      /|         / |
 *    5/_|________/6 |
 *    |  |        |  |
 *    |  |        |  |
 *    | 0|________|__|3
 *    | /         | /
 *  1 |/__________|/ 2      ->x
 *
 *  /
 * v Z
 */
class Hexahedral : public Element
{
public:
    /**
     * @brief Конструктор
     * @param node_0 Номер нулевого узла
     * @param node_1 Номер первого узла
     * @param node_2 Номер второго узла
     * @param node_3 Номер третьего узла
     * @param node_4 Номер четвертого узла
     * @param node_5 Номер пятого узла
     * @param node_6 Номер шестого узла
     * @param node_7 Номер седьмого узла
     */
    Hexahedral(const UInteger &node_0, const UInteger &node_1, const UInteger &node_2, const UInteger &node_3, const UInteger &node_4, const UInteger &node_5, const UInteger &node_6, const UInteger &node_7);
    /**
     * @brief Конструктор копирования
     * @param hexahedral Объект-шестигранник для копирования
     */
    Hexahedral(const Hexahedral &hexahedral);
    /**
     * @brief Оператор доступа к вершине по индексу
     * @param i Индекс вершины шестигранника
     * @return Номер узла, на который ссылается вершина
     */
    UInteger &operator [] (int i);
    /**
     * @brief Оператор присвоения
     * @param hexahedral Объект-шестигранник для копирования
     * @return Копию объекта hexahedral
     */
    Hexahedral &operator =(const Hexahedral &hexahedral);
    /**
     * @brief Количество узлов
     * @return Количество узлов (вершин) в элементе (8)
     */
    virtual int verticesCount () const;
    /**
     * @brief Номер узла, на который указывает вершина элемента
     * @param i Номер вершины
     * @return Номер узла, соответствующего вершине с номером i
     */
    virtual UInteger vertexNode (int i) const;
    /**
     * @brief Количество граней, соответвующих элементу (необходимо для визуализации)
     * @return Количество граней (6)
     */
    virtual int facesCount() const;
    /**
     * @brief face Грань элемента
     * @param i Номер грани
     * @return Грань элемента, номер которойс равен i
     */
    virtual UIntegerVector face(const int &i) const;
    /**
     * @brief Проверка принадлежности узла элементу
     * @param node Номер (хеш-код) узла
     * @return true, если узел принадлежит элементу, false, иначе
     */
    virtual bool in(const UInteger &node) const;
    /**
     * @brief Возвращает индекс узла, номер (хеш-код) котороый равен заданному
     * @param node Номер (хеш-код) узла
     * @return Индекс (значение от 0 до verticesCount() - 1) если такой номер найден, -1 - иначе
     */
    virtual int index(const UInteger &node) const;
private:
    UInteger vertex_[8]; //!< Вершины шестигранника
};
}


#endif // HEXAHEDRAL_H
