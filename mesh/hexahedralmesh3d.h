/**
  * @author Сергей Чопоров
  * @date 23/06/2014
  * @version 1.0.1
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
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
     * @brief Метод создает равномерную сетку шестигранников для прямоугольной призмы
     * @param xCount Количество точек дискретизации по оси ординат
     * @param yCount Количество точек дискретизации по оси абсцисс
     * @param zCount Количество точек дискретизации по оси аппликат
     * @param xMin Абсцисса минимального угла
     * @param yMin Ордината минимального угла
     * @param zMin Аппликата минимального угла
     * @param width Ширина
     * @param height Высота
     * @param depth Глубина
     */
    void prismDomain(const UInteger &xCount, const UInteger &yCount, const UInteger &zCount,
                     const double &xMin, const double &yMin, const double &zMin,
                     const double &width, const double &height, const double &depth);
    /**
     * @brief Метод строит сетку как тело вращения
     * @param baseMesh Указатель на базовую двумерную сетку четырехугольных элементов
     * @param xDelta Смещение вдоль оси ординат
     * @param yDelta Смещение вдоль оси абсцисс
     * @param lCount Количество слоев элементов
     * @param x_axes Ось вращения (x_axes == true, то вращение вокруг оси ординат; иначе - абсцисс)
     * @param withLayersInfo Учеть номера слоев
     */
    void rotateBaseMesh(QuadrilateralMesh2D *baseMesh, const double& xDelta, const double& yDelta, const int& lCount, bool x_axes, bool withLayersInfo = true);
    /**
     * @brief Метод строит сетку как тело вращения на заданный угол
     * @param baseMesh Указатель на базовую двумерную сетку четырехугольных элементов
     * @param xDelta Смещение вдоль оси ординат
     * @param yDelta Смещение вдоль оси абсцисс
     * @param angle Угол поворота профиля
     * @param lCount Количество слоев элементов
     * @param x_axes Ось вращения (x_axes == true, то вращение вокруг оси ординат; иначе - абсцисс)
     * @param withLayersInfo Учеть номера слоев
     */
    void rotateBaseMesh(QuadrilateralMesh2D *baseMesh, const double& xDelta, const double& yDelta, const double& angle, const int& lCount, bool x_axes, bool withLayersInfo = true);
    /**
     * @brief Конструктор копирования
     * @param mesh Экземпляр объекта для копирования
     */
    HexahedralMesh3D(const HexahedralMesh3D &mesh);
    /**
     * @brief Конструктор создает копию объекта, переданого по ссылке
     * @param mesh Ссылка на объект для копирования
     */
    HexahedralMesh3D(const HexahedralMesh3D *mesh);
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
    /**
     * @brief Вычислить площадь грани (грань - четырехугольник)
     * @param face Список номеров узлов, определяющих грань
     * @return Площадь грани
     */
    double faceArea(const UIntegerVector &face) const;
    /**
     * @brief Вычислить площадь поверхности дискретной модели
     * @return Площадь поверхности дискретной модели
     */
    virtual double surfaceArea() const;
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
    /**
     * @brief Добавить элемент к сетке
     * @param Массив ссылок (номеров) на узлы
     */
    void addElement(const std::vector<UInteger> &nodes_ref);
    /**
     * @brief Метод очищает информацию об елементах
     */
    virtual void clearElements();
    /**
     * @brief Генерация сетки движением вдоль оси z
     * @param baseMesh Указатель на базовую сетку (профиль)
     * @param z0 Ближння координата профиля
     * @param z1 Дальняя координата профиля
     * @param phi0 Начальный поворот профиля
     * @param phi1 Конечный поворот профиля
     * @param k0 Начальный масштаб профиля
     * @param k1 Конечный масштаб профиля
     * @param zLayersCount Количество слоев вдоль оси z
     */
    void sweepBaseMesh(QuadrilateralMesh2D *baseMesh, const double &z0, const double &z1, const double &phi0, const double &phi1, const double &k0, const double &k1, const int &zLayersCount);
//    UInteger toArray(UInteger i, UInteger j, UInteger k, UInteger yCount, UInteger zCount) { return i * yCount * zCount + j * zCount + k; }
private:
    std::vector<Hexahedral> element_; //!< Массив шестигранных элементов
};
}


#endif // HEXAHEDRALMESH3D_H
