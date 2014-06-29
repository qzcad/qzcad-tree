/**
  * @author Сергей Чопоров
  * @date 15/10/2012
  * @version 1.0.5
  */
#ifndef ELEMENT_H
#define ELEMENT_H
#include "integer.h"
#include "integervector.h"
namespace msh
{
/**
 * @brief Класс Element - абстаркция элемента сетки
 */
class Element
{
public:
    virtual ~Element() {}
    /**
     * @brief Количество узлов
     * @return Количество узлов (вершин) в элементе
     */
    virtual int verticesCount () const = 0;
    /**
     * @brief Номер узла, на который указывает вершина элемента
     * @param i Номер вершины
     * @return Номер узла, соответствующего вершине с номером i
     */
    virtual UInteger vertexNode (int i) const = 0;
    /**
     * @brief Количество граней, соответвующих элементу (необходимо для визуализации)
     * @return Количество граней
     */
    virtual int facesCount() const = 0;
    /**
     * @brief face Грань элемента
     * @param i Номер грани
     * @return Грань элемента, номер которойс равен i
     */
    virtual UIntegerVector face(const int &i) const = 0;
};
}

#endif // ELEMENT_H
