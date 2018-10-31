#ifndef TETRAHEDRON_H
#define TETRAHEDRON_H

#include "element.h"

namespace msh {
/**
 * @brief Класс Tetrahedron - четырехугольный тетраэдр
 * ^
 * | Y
 *    3____0
 *    |%  /\
 *    | %/  \
 *    | / %  \
 *    |/    % \
 *  1 |_______%\ 2   ->x
 *
 *  /
 * v Z
 */
class Tetrahedron : public Element
{
public:
    /**
     * @brief Конструктор
     * @param node_0 Номер нулевого узла
     * @param node_1 Номер первого узла
     * @param node_2 Номер второго узла
     * @param node_3 Номер третьего узла
     */
    Tetrahedron(const UInteger &node_0, const UInteger &node_1, const UInteger &node_2, const UInteger &node_3);
    /**
     * @brief Конструктор копирования
     * @param tetrahedron Объект тераэдр для копирования
     */
    Tetrahedron(const Tetrahedron &tetrahedron);
    /**
     * @brief Оператор доступа к вершине по индексу
     * @param i Индекс вершины тетраэдра
     * @return Номер узла, на который ссылается вершина
     */
    UInteger &operator [] (int i);
    /**
     * @brief Константный оператор доступа к вершине по индексу
     * @param i Индекс вершины тетраэдра
     * @return Номер узла, на который ссылается вершина
     */
    const UInteger &operator [](int i) const;
    /**
     * @brief Оператор присвоения
     * @param hexahedral Объект-тетраэдр для копирования
     * @return Копию объекта hexahedral
     */
    Tetrahedron &operator =(const Tetrahedron &tetrahedron);
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
     * @return Количество гранейn (6)
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
    UInteger vertex_[4]; //!< Вершины тетраэдра
};
}

#endif // TETRAHEDRON_H
