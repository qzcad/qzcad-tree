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
class TetrahedralMesh3D;
}
#include "trianglemesh3d.h"
#include "hexahedralmesh3d.h"

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
     * @brief Конструктор создает копию объекта, переданого по указателю
     * @param mesh Указатель на объект для копирования
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
    void addElement(const Tetrahedron &t);
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
    void addElement(const UInteger &i0, const UInteger &i1, const UInteger &i2, const UInteger &i3, const UInteger &i4, const UInteger &i5);
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
    void sweepBaseMesh(TriangleMesh2D *baseMesh, const double &z0, const double &z1, const double &phi0, const double &phi1, const double &k0, const double &k1, const int &zLayersCount);
    /**
     * @brief Конвертировать сетку шестигранников (используя схему без вставки дополнительных узлов)
     * @param Указатель на сетку шестигранников
     */
    void convertHexahedralMesh(const HexahedralMesh3D *mesh);
    /**
     * @brief Конвертировать сетку шестигранников (используя схему Face-Centred: 24 тетраэдра на 1 шестигранник)
     * @param Указатель на сетку шестигранников
     */
    void convertHexahedralMeshFC(const HexahedralMesh3D *mesh);
    /**
     * @brief Напечать статистику дискретной модели
     */
    void printStats() const;
    /**
     * @brief Триангуляция Делоне контура функциональной модели
     * @param mesh Объект-конутр
     * @param func Указатель на функцию области. Если он равен nullptr, то используется сfunction.
     * @see TriangleMesh3D::cfunction
     */
    void delaunay(const TriangleMesh3D &mesh, std::function<double(double, double, double)> func);
    /**
     * @brief Метод построения стеки с использованием фоновых тетраэдров
     * @param mesh Сетка шестигранников
     * @param func Функция области
     * @param level Уровень нуля
     * @param smooth Количество итераций сглаживания
     * @param optimize Количество итераций оптимизации
     */
    void backgroundGrid(const TetrahedralMesh3D *mesh, std::function<double(double, double, double)> func, double level = 0.0, int smooth = 0, int optimize = 0);
    void localFubctionalOptimization(int maxiter=2);
    /**
     * @brief Вычислить соотношение длин сторон (минимальной к максимальной)
     * @param elnum Номер элемента
     * @return Соотношение длин сторон (минимальной к максимальной)
     */
    virtual double lengthAspect(const UInteger &elnum) const;
    /**
     * @brief Вычислить значение минимального угла в элементе
     * @param elnum Номер элемента
     * @return Минимальный угол элемента (радианы)
     */
    virtual double minAngle(const UInteger &elnum) const;
    /**
     * @brief Вычислить значение максимального угла в элементе
     * @param elnum Номер элемента
     * @return Минимальный угол элемента (радианы)
     */
    virtual double maxAngle(const UInteger &elnum) const;
    void subdivide(std::list<UInteger> eNumbers, std::function<double(double, double, double)> func);
private:
    std::vector<Tetrahedron> element_; //!< Массив шестигранных элементов
};
}

#endif // TETRAHEDRALMESH3D_H
