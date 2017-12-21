/**
  * @author Сергей Чопоров
  * @date 23/01/2014
  * @version 1.0.1
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef SEGMENT_H
#define SEGMENT_H

#include "element.h"

namespace msh
{
/**
 * @brief Класс Segment - сегмент контура
 */
class Segment : public Element
{
public:
    /**
     * @brief Segment Конструктор
     * @param left Левая (первая) вершина
     * @param right Правая (вторая) вершина
     */
    Segment(const long &left, const long &right);
    /**
     * @brief Segment Конструктор копирования
     * @param segment Сегмент-объект для копирования
     */
    Segment(const Segment &segment);
    virtual ~Segment() {}
    /**
     * @brief verticesCount Количество вершин в сегменте
     * @return 2
     */
    int verticesCount () const;
    /**
     * @brief vertexNode Номер узла, на который указывает вершина сегмента
     * @param i Номер вершины
     * @return Номер узла, соответствующего вершине с номером i
     */
    UInteger vertexNode (int i) const;
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
     * @brief Проверка на совпадение пары соединенных узлов (без учета направления)
     * @param s Отрезок, с которым необходимо сравниться
     * @return ,если пара соединенных узлов совпадает
     */
    bool isSame(const Segment &s);
    /**
     * @brief Проверка принадлежности узла элементу
     * @param node Номер (хеш-код) узла
     * @return true, если узел принадлежит элементу, false, иначе
     */
    virtual bool in(const UInteger &node) const;
private:
    UInteger leftVertex_; //!< левая вершина
    UInteger rightVertex_; //!< правая вершина
};
}

#endif // SEGMENT_H
