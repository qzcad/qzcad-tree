/**
  * @author Сергей Чопоров
  * @date 05/08/2015
  * @version 1.0.0
  * @copyright Copyright 2015 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef QUADRILATERALMESH3D_H
#define QUADRILATERALMESH3D_H

#include <functional>

#include "mesh3d.h"
#include "quadrilateral.h"
#include "hexahedralmesh3d.h"

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
     * @brief Конструктор создает копию объекта, переданого по указателю
     * @param mesh Указатель на объект для копирования
     */
    QuadrilateralMesh3D(const QuadrilateralMesh3D *mesh);
    /**
     * @brief Метод для создания равномерной сетки в цилиндрических координатах
     * @param rCount Количество элементов вдоль радиуса
     * @param lCount Количество элементов вдоль образующей
     * @param radius Радиус
     * @param length Длина образующей
     */
    void cylinderDomain(const UInteger &rCount, const UInteger &lCount, const double &radius, const double &length);
    void cylinderDomain(const UInteger &rCount, const UInteger &lCount, std::function<double (double)> radius, const double &length);
    /**
     * @brief Метод для создания равномерной сетки в конических координатах
     * @param rCount Количество элементов вдоль радиуса
     * @param lCount Количество элементов вдоль образующей
     * @param bottom_radius Нижний радиус
     * @param top_radius Верхний радиус
     * @param length Длина образующей
     */
    void coneDomain(const UInteger &rCount, const UInteger &lCount, const double &bottom_radius, const double &top_radius, const double &length);
    /**
     * @brief Метод построения стеки с использованием фоновых шестигранников
     * @param mesh Сетка шестигранников
     * @param func Функция области
     * @param level Уровень нуля
     * @param smooth Количество итераций сглаживания
     * @param optimize Количество итераций оптимизации
     */
    void backgroundGrid(const HexahedralMesh3D *mesh, std::function<double(double, double, double)> func, double level = 0.0, int smooth = 0, int optimize = 0);
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
     * @brief Добавить элемент к сетке
     * @param Массив ссылок (номеров) на узлы
     */
    void addElement(const std::vector<UInteger> &nodes_ref);
    /**
     * @brief area Подсчитать площадь элемента
     * @param number Номер элемента
     * @return Площадь элемента
     */
    virtual double area(const UInteger &number) const;
    /**
     * @brief Операция добавления сетки к существующей на основе элементарного объединения
     * @param mesh Указатель на сетку для добавления
     */
    void add(const QuadrilateralMesh3D *mesh);
    /**
     * @brief Операция перемещения сетки на заданный радиус вектор
     * @param x Абсцисса радиус вектора перемещения
     * @param y Ордината радиус вектора перемещения
     * @param y Аппликата радиус вектора перемещения
     */
    void translate(const double &x, const double &y, const double &z);
    /**
     * @brief Метод очищает информацию об елементах
     */
    virtual void clearElements();
    /**
     * @brief Процедура сглаживания Лапласа
     * @param func Функция области
     * @param level Линия уровня
     * @param iter_num Количесво итераций
     */
    void laplacianSmoothing(std::function<double(double, double, double)> func, double level = 0, int iter_num = 1);
    /**
     * @brief Процедура сглаживания на основе минимизации функционала расстояния-длины
     * @param func Функция области
     * @param level Линия уровня
     * @param iter_num Количесво итераций
     */
    void distlenSmoothing(std::function<double(double, double, double)> func, double level = 0, int iter_num = 1);
    /**
     * @brief Вычислить значение минимального угла в элементе
     * @param elNum Номер элемента
     * @return Минимальный угол элемента (радианы)
     */
    double minAngle(const UInteger &elNum) const;
    /**
     * @brief Вычислить значение максимального угла в элементе
     * @param elNum Номер элемента
     * @return Минимальный угол элемента (радианы)
     */
    double maxAngle(const UInteger &elNum) const;
    double maxAngle(const Point3D &A, const Point3D &B, const Point3D &C, const Point3D &D) const;
    void flip();
    void refineTopology(std::function<double(double, double, double)> func, double level = 0);
protected:
    std::vector<Quadrilateral> element_; //!< Массив элементов
};
}

#endif // QUADRILATERALMESH3D_H
