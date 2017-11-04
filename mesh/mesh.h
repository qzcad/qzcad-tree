  /**
  * @author Сергей Чопоров
  * @date 15/10/2012
  * @version 1.0.6
  * @copyright Copyright 2012 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef MESH_H
#define MESH_H

#include <vector>
#include <string>

#include "pointpointer.h"
#include "elementpointer.h"
#include "nodetype.h"
#include "nameddoublevector.h"
#include "adjacentset.h"
#include "point3d.h"

namespace msh
{
/**
 * @brief Класс Mesh - абстракция сетки элементов
 * В классе Mesh не определены тип сетки и форма элемента
 * @see Element, Point, UInteger, Floating, NodeType
 */
class Mesh
{

public:
    /**
     * @brief Конструктор, обеспечивающий копирование всей дополнительной информации
     * @param mesh
     */
    Mesh(const Mesh *mesh);
    /**
     * @brief Количество узлов
     * @return Количество узлов в сетке
     */
    virtual UInteger nodesCount() const = 0;
    /**
     * @brief Координаты узла сетки
     * @param number Номер узла сетки
     * @return Указатель на координаты узла сетки с номером number
     */
    virtual PointPointer node(const UInteger &number) const = 0;
    /**
     * @brief Тип узла сетки
     * @param number Номер узла сетки
     * @return Тип узла сетки с номером number
     */
    virtual NodeType nodeType(const UInteger &number) const = 0;
    /**
     * @brief Количество элементов
     * @return Количество элементов в сетке
     */
    virtual UInteger elementsCount() const = 0;
    /**
     * @brief Элемет сетки
     * @param number Номер элемента сетки
     * @return Указатель на элемент сетки
     */
    virtual ElementPointer element(const UInteger &number) const = 0;
    /**
     * @brief Минимальное значение ординаты точек сетки
     * @return Минимальное значение ординаты точек сетки
     */
    virtual double xMin() const = 0;
    /**
     * @brief Максимальное значение ординаты точек сетки
     * @return Максимальное значение ординаты точек сетки
     */
    virtual double xMax() const = 0;
    /**
     * @brief Минимальное значение абсциссы точек сетки
     * @return Минимальное значение абсциссы точек сетки
     */
    virtual double yMin() const = 0;
    /**
     * @brief Максимальное значение абсциссы точек сетки
     * @return Максимальное значение абсциссы точек сетки
     */
    virtual double yMax() const = 0;
    /**
     * @brief Минимальное значение аппликаты точек сетки
     * @return Минимальное значение аппликаты точек сетки
     */
    virtual double zMin() const = 0;
    /**
     * @brief Максимальное значение аппликаты точек сетки
     * @return Максимальное значение аппликаты точек сетки
     */
    virtual double zMax() const = 0;
    /**
     * @brief Размерность пространства
     * @return Размерность пространства, в котором определена сетка
     */
    virtual int dimesion() const = 0;
    /**
     * @brief Определить принадлежность элемента границе
     * @param element Константная ссылка на элемент
     * @return true - граничный элемент; false - внутренний
     */
    virtual bool isBorderElement(ElementPointer element) const;
    /**
     * @brief Определить принадлежность грани элемента границе
     * @param face Масив, представляющий грань
     * @return true - граничная грань; false - внутреняя
     */
    virtual bool isBorderFace(const UIntegerVector &face) const;
    /**
     * @brief Вычислить нормаль к грани
     * @param face Масив, представляющий грань
     * @return Координаты вектора-нормали
     */
    virtual Point3D normal(const UIntegerVector &face) const = 0;
    /**
     * @brief Обновить параметры области определения сетки (xMin, xMax, yMin, yMax, zMin, zMax)
     */
    virtual void updateDomain() = 0;
    /**
     * @brief Метод возвращает множество номеров смежных в узле элементов
     * @param nodeNumber Номер узла
     * @return Множество смежных элементов (их номера)
     */
    virtual AdjacentSet adjacent(const UInteger &nodeNumber) const = 0;
    /**
     * @brief Количество элементов, смежных в узле с заданным номером
     * @param nodeNumber Номер узла
     * @return Количество элементов, смежных в узле
     */
    virtual UInteger adjacentCount(const UInteger &nodeNumber) const = 0;
    /**
     * @brief Добавить узел в сетку
     * @param point Указатель на координаты
     * @param type Тип узла
     * @return Номер вставленного узла в сетке
     */
    virtual UInteger pushNode(PointPointer point, const NodeType &type) = 0;
    /**
     * @brief Добавить элемент к сетке
     * @param Массив ссылок (номеров) на узлы
     */
    virtual void addElement(const std::vector<UInteger> &nodes_ref) = 0;
    /**
     * @brief Номер слоя элемента с заданным номером
     * @param number Номер элемента
     * @return Номер слоя для элемента с заданным номером
     */
    virtual int layer(const UInteger &number) const;
    /**
     * @brief Установить номер слоя элемента
     * @param number Номер элемента
     * @param l Номер слоя/
     */
    virtual void setLayer(const UInteger &number, const int &l);
    /**
     * @brief pushLayer Добавить значение с номером слоя в массив значений
     * @param l Номер слоя
     */
    virtual void pushLayer(const int &l);
    /**
     * @brief Очистить информацию о номрах слоя
     */
    virtual void clearLayers();
    /**
     * @brief sizeOfLayers
     */
    virtual UInteger sizeOfLayers() const;
    /**
     * @brief Виртуальный деструктор
     */
    virtual ~Mesh(){}
    /**
     * @brief Получить значение текущего уровня нуля (точность вычислительных операций)
     * @return Значение текущего уровня нуля
     */
    static double epsilon();
    /**
     * @brief Установить значение текущего уровня нуля (точность вычислительных операций)
     * @param epsilon Новое значение текущего уровня нуля
     */
    static void setEpsilon(double epsilon);
    /**
     * @brief Метод для получения количества массивов значений с дополнительной числовой инфориацией
     * @return Количество массивов значений с дополнительной числовой инфориацией
     */
    UInteger dataVectorsCount() const;
    /**
     * @brief Добавить именованный вектор значений
     * @param vec Вектор значений
     */
    void addDataVector(const NamedDoubleVector &vec);
    /**
     * @brief Добавить именованный вектор значений
     * @param name Название вектора
     * @param values Значения
     */
    void addDataVector(const std::string &name, const std::vector<double> &values);
    /**
     * @brief Метод для получения массива значений с заданным номером
     * @param i Номер массива
     * @return Массив значений с заданным номером
     */
    const NamedDoubleVector &data(const UInteger &i) const;
    /**
     * @brief Метод очищает массивы значений
     */
    virtual void clearDataVectors();
    /**
     * @brief Метод очищает информацию об узлах сетки
     */
    virtual void clearNodes() = 0;
    /**
     * @brief Метод очищает информацию об елементах
     */
    virtual void clearElements() = 0;
    /**
     * @brief Метод очищает дискретную модель
     */
    virtual void clear();
    /**
     * @brief Напечатать в стандартный вывод экстремальные значения векторов данных
     */
    void printDataExtremums();
    /**
     * @brief Определить степень грани
     * @param face Масив, представляющий грань
     * @return Степень грани (число смежных элементов)
     */
    virtual int facePower(const UIntegerVector &face) const;

protected:
    std::vector<int> layer_; //!< Массив с номером слоя для каждого элемента
    static double epsilon_; //!< Точность вычислительных операций, по умолчантю 1.0E-10
    std::vector<NamedDoubleVector> data_; //!< Дополнительная числовая информация, определенная на сетке
};
}
#endif // MESH_H
