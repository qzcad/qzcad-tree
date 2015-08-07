/**
  * @author Сергей Чопоров
  * @date 05/08/2015
  * @version 1.0.0
  * @copyright Copyright 2015 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef QUADRILATERALMESH3D_H
#define QUADRILATERALMESH3D_H

#include "mesh3d.h"
#include "quadrilateral.h"

namespace msh {
/**
 * @brief Класс QuadrilateralMesh3D - абстракция сетки четырехугольных элементов в трехмерном пространстве (поверхностная сетка).
 */
class QuadrilateralMesh3D : public Mesh3D
{
public:
    /**
     * @brief Конструктор по умолчанию
     */
    QuadrilateralMesh3D();
    /**
     * @brief Конструктор копирования
     * @param mesh Экземпляр объекта для копирования
     */
    QuadrilateralMesh3D(const QuadrilateralMesh3D &mesh);
    /**
     * @brief Конструктор создает копию объекта, переданого по ссылке
     * @param mesh Ссылка на объект для копирования
     */
    QuadrilateralMesh3D(const QuadrilateralMesh3D *mesh);
    /**
     * @brief Конструктор для создания равномерной сетки в цилиндрических координатах
     * @param rCount Количество элементов вдоль радиуса
     * @param lCount Количество элементов вдоль образующей
     * @param radius Радиус
     * @param length Длина образующей
     */
    QuadrilateralMesh3D(const UInteger &rCount, const UInteger &lCount, const double &radius, const double &length);
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
     * @brief Определить принадлежность элемента границе
     * @param number Номер элемента
     * @return true - граничный элемент; false - внутренний
     */
    virtual bool isBorderElement(const UInteger &number) const;
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
     * Этот метод не является потокобезопасным.
     */
    void addElement(const UInteger &node0, const UInteger &node1, const UInteger &node2, const UInteger &node3);
    /**
     * @brief Добавить элемент к сетке
     * @param quad Ссылка на элемент, который необходимо добавить к сетке
     * Этот метод не является потокобезопасным.
     */
    void addElement(const Quadrilateral &quad);
    /**
     * @brief area Подсчитать площадь элемента
     * @param number Номер элемента
     * @return Площадь элемента
     */
    virtual double area(const UInteger &number) const;
protected:
    std::vector<Quadrilateral> element_; //!< Массив элементов
};
}

#endif // QUADRILATERALMESH3D_H