/**
  * @author Сергей Чопоров
  * @date 16/06/2015
  * @version 1.0.0
  * @copyright Copyright 2015 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  * */
#ifndef GAUSSQUADRATURE_H
#define GAUSSQUADRATURE_H

#include "doublevector.h"

using namespace mtx;

/**
 * @brief Класс для генерации квадратур Гаусса
 */
class GaussQuadrature
{
public:
    /**
     * @brief Конструктр по умолчанию
     */
    GaussQuadrature();
    /**
     * @brief Метод для генерации квадратур для интегрирования на отрезке [-1; 1]
     * @param count Количество точек в квадратуре (от 1 до 5)
     * @param point Массив для записи координат на отрезке
     * @param weight Массив для записи весовых коэффициентов соответствующих координат
     */
    void quadrature(int count, DoubleVector &point, DoubleVector &weight);
    /**
     * @brief Метод для генерации квадратур для интегрирования на отрезке [0; 1] для интегрирования на треугольнике
     * @param count Количество точек в квадратуре (1, 3, 4)
     * @param xi Координата первого направления
     * @param eta Координата второго направления
     * @param weight Массив для записи весовых коэффициентов соответствующих координат
     */
    void quadrature(int count, DoubleVector &xi, DoubleVector &eta, DoubleVector &weight);
};

#endif // GAUSSQUADRATURE_H
