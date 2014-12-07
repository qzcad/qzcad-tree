/**
  * @author Сергей Чопоров
  * @date 15/10/2012
  * @version 1.0.3
  * @copyright Copyright 2012 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef POINT_H
#define POINT_H

namespace msh
{
/**
 * @brief Класс Point - абстракция точки в пространстве
 * Одномерные и двумерные точки можно рассматривать как проекции или
 * как трехмерные точки с соответствующими координатами равными нулю
 * Класс описывает интерфейс доступа к координатам точки
 */
class Point
{
public:
    virtual ~Point() {}
    /**
     * @brief Размерность пространства
     * @return размерность пространства, в ктором определена точка
     */
    virtual int dimension() const = 0;
    /**
     * @brief Ордината
     * @return Ординату точки или ее проекции
     */
    virtual double x() const = 0;
    /**
     * @brief Абсцисса
     * @return Абсциссу точки или ее проекции
     */
    virtual double y() const = 0;
    /**
     * @brief Аппликата
     * @return Аппликату точки или ее проекции
     */
    virtual double z() const = 0;
    /**
     * @brief Длина вектора
     * @return Длину вектора, определенного началом координат и координатами точки
     */
    virtual double length() const = 0;
};
}

#endif // POINT_H
