/**
  * @author Сергей Чопоров
  * @date 23/06/2014
  * @version 1.0.1
  */
#ifndef HEXAHEDRALMESH3D_H
#define HEXAHEDRALMESH3D_H

#include <vector>

#include "mesh3d.h"
#include "hexahedral.h"
#include "quadrilateralmesh2d.h"

namespace msh {
/**
 * @brief Класс HexahedralMesh3D - сетка шестигранных элементов в трехмерном пространстве
 * @see Mesh3D, Hexahedral
 */
class HexahedralMesh3D : public Mesh3D
{
public:
    HexahedralMesh3D();
    /**
     * @brief Конструктор, который создает равномерную сетку шестигранников для прямоугольной призмы
     * @param xCount Количество точек дискретизации по оси ординат
     * @param yCount Количество точек дискретизации по оси абсцисс
     * @param zCount Количество точек дискретизации по оси аппликат
     * @param xMin Ордината минимального угла
     * @param yMin Абсцисса минимального угла
     * @param zMin Аппликата минимального угла
     * @param width Ширина
     * @param height Высота
     * @param depth Глубина
     */
    HexahedralMesh3D(const UInteger &xCount, const UInteger &yCount, const UInteger &zCount,
                     const Floating &xMin, const Floating &yMin, const Floating &zMin,
                     const Floating &width, const Floating &height, const Floating &depth);
    /**
     * @brief Конструктор, который строит сетку как тело вращения
     * @param baseMesh Указатель на базовую двумерную сетку четырехугольных элементов
     * @param xDelta Смещение вдоль оси ординат
     * @param yDelta Смещение вдоль оси абсцисс
     * @param lCount Количество слоев элементов
     * @param x_axes Ось вращения (x_axes == true, то вращение вокруг оси ординат; иначе - абсцисс)
     * @param withElementValue Учеть значения, опредедленные на элементе
     */
    HexahedralMesh3D(QuadrilateralMesh2D *baseMesh, const Floating& xDelta, const Floating& yDelta, const int& lCount, bool x_axes, bool withElementValue = true);
    /**
     * @brief HexahedralMesh3D Конструктор, который строит сетку как тело вращения на заданный угол
     * @param baseMesh Указатель на базовую двумерную сетку четырехугольных элементов
     * @param xDelta Смещение вдоль оси ординат
     * @param yDelta Смещение вдоль оси абсцисс
     * @param angle Угол поворота профиля
     * @param lCount Количество слоев элементов
     * @param x_axes Ось вращения (x_axes == true, то вращение вокруг оси ординат; иначе - абсцисс)
     * @param withElementValue Учеть значения, опредедленные на элементе
     */
    HexahedralMesh3D(QuadrilateralMesh2D *baseMesh, const Floating& xDelta, const Floating& yDelta, const Floating& angle, const int& lCount, bool x_axes, bool withElementValue = true);
    /**
     * @brief Количество элементов
     * @return Количество элементов в сетке
     */
    virtual UInteger elementsCount() const;
    /**
     * @brief Элемет сетки
     * @param number Номер элемента сетки
     * @return Указатель на копию элемента сетки
     */
    virtual ElementPointer element(const UInteger &number) const;
    virtual bool isBorderElement(const UInteger &number) const;
    /**
     * @brief Значение некоторой функции, определенной на элементе
     * @param number Номер элемента
     * @return Значение, соответствующее элементу
     */
    virtual Floating elementValue(const UInteger &number) const;
protected:
    /**
     * @brief Добавить элемент к  сетке
     * @param node0 Номер узла в вершине 0
     * @param node1 Номер узла в вершине 1
     * @param node2 Номер узла в вершине 2
     * @param node3 Номер узла в вершине 3
     * @param node4 Номер узла в вершине 4
     * @param node5 Номер узла в вершине 5
     * @param node6 Номер узла в вершине 6
     * @param node7 Номер узла в вершине 7
     */
    void addElement(const UInteger &node0, const UInteger &node1, const UInteger &node2, const UInteger &node3,
                    const UInteger &node4, const UInteger &node5, const UInteger &node6, const UInteger &node7);
//    UInteger toArray(UInteger i, UInteger j, UInteger k, UInteger yCount, UInteger zCount) { return i * yCount * zCount + j * zCount + k; }
private:
    std::vector<Hexahedral> element_; //!< Массив шестигранных элементов
    std::vector<Floating> elementValue_; //!< Значение на элементе
};
}


#endif // HEXAHEDRALMESH3D_H
