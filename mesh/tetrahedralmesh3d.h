/**
  * @author Сергей Чопоров
  * @date 14/03/2017
  * @version 1.0.0
  * @copyright Copyright 2017 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef TETRAHEDRALMESH3D_H
#define TETRAHEDRALMESH3D_H

#include <vector>

#include "mesh3d.h"
#include "tetrahedron.h"
#include "trianglemesh2d.h"

namespace msh {
/**
 * @brief Класс TetrahedralMesh3D - сетка тетраэдров в трехмерном пространстве
 * @see Mesh3D, Tetrahedron
 */
class TetrahedralMesh3D : public Mesh3D
{
public:
    TetrahedralMesh3D();
    /**
     * @brief Конструктор копирования
     * @param mesh Экземпляр объекта для копирования
     */
    TetrahedralMesh3D(const TetrahedralMesh3D &mesh);
    /**
     * @brief Конструктор создает копию объекта, переданого по ссылке
     * @param mesh Ссылка на объект для копирования
     */
    TetrahedralMesh3D(const TetrahedralMesh3D *mesh);
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
     */
    void addElement(const UInteger &node0, const UInteger &node1, const UInteger &node2, const UInteger &node3);
    /**
     * @brief Добавить элемент к сетке
     * @param Массив ссылок (номеров) на узлы
     */
    void addElement(const std::vector<UInteger> &nodes_ref);
    /**
     * @brief Метод очищает информацию об елементах
     */
    virtual void clearElements();
    void sweepBaseMesh(TriangleMesh2D *baseMesh, const double &z0, const double &z1, const double &phi0, const double &phi1, const double &k0, const double &k1, const int &zLayersCount);
private:
    std::vector<Tetrahedron> element_; //!< Массив шестигранных элементов
};
}

#endif // TETRAHEDRALMESH3D_H
