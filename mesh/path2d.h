/**
  * @author Сергей Чопоров
  * @date 21/01/2014
  * @version 1.0.1
  */
#ifndef PATH2D_H
#define PATH2D_H

#include <vector>

#include "mesh.h"
#include "point2d.h"

namespace msh
{
/**
 * @brief Класс Path2D - контур в двумерном пространстве
 */
class Path2D : public Mesh
{
public:
    /**
     * @brief Path2D Конструктор по умолчанию
     */
    Path2D();
    /**
     * @brief Path2D Конструктор копирования
     * @param path2d Объект для копирования
     */
    Path2D(const Path2D &path2d);

    double xMin() const;
    void setXMin(const double &xMin);

    double xMax() const;
    void setXMax(const double &xMax);

    double yMin() const;
    void setYMin(const double &yMin);

    double yMax() const;
    void setYMax(const double &yMax);
    /**
     * @brief nodesCount Количество узлов в контуре
     * @return Количество узлов, образующих контур
     */
    UInteger nodesCount() const;
    /**
     * @brief node Координаты узла
     * @param number Номер узла
     * @return Координаты узла с номером number
     */
    PointPointer node(UInteger number) const;
    /**
     * @brief elementsCount Количество элементов
     * @return Количество отрезков (граней), которые образуют контур
     */
    UInteger elementsCount() const;
    /**
     * @brief element Элемент сетки
     * @param number Номер элемента сетки
     * @return Указатель на элемент - отрезок контура
     */
    Element *element(UInteger number) const;
    /**
     * @brief zMin Минимальное значение аппликаты точек сетки
     * @return Минимальное значение аппликаты точек сетки
     */
    double zMin() const;
    /**
     * @brief zMax Максимальное значение аппликаты точек сетки
     * @return Максимальное значение аппликаты точек сетки
     */
    double zMax() const;
    /**
     * @brief dimesion Размерность пространства
     * @return 2
     */
    int dimesion() const;
    /**
     * @brief addNode Добавить узел
     * @param point Координаты узла
     */
    void addNode(const Point2D &point);
    /**
     * @brief removeNode Удалить узел с заданным номером
     * @param number Номеро узла, который необходимо удалить из контура
     */
    void removeNode(UInteger number);
    /**
     * @brief insertNode Вставить узел перед узлом с заданным номером
     * @param number Номер узла, перед которым необходимо произвести вставку
     */
    void insertNode(UInteger number, const Point2D &point);

private:
    std::vector<Point2D> node_; //!< Список узлов
    double xMin_; //!< Минимальное значение ординаты
    double xMax_; //!< Максимальное значение ординаты
    double yMin_; //!< Минимальное значение абсциссы
    double yMax_; //!< Максимальное значение абсциссы
};
}

#endif // PATH2D_H
