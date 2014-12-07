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

private:
    UInteger leftVertex_; //!< левая вершина
    UInteger rightVertex_; //!< правая вершина
};
}

#endif // SEGMENT_H
