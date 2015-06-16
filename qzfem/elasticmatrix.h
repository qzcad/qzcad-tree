/**
  * @author Сергей Чопоров
  * @date 15/06/2015
  * @version 1.0.0
  * @copyright Copyright 2015 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  * */
#ifndef ELASTICMATRIX_H
#define ELASTICMATRIX_H

#include "doublematrix.h"

using namespace mtx;
/**
 * @brief Класс для построения матрицы упругости
 */
class ElasticMatrix
{
public:
    /**
     * @brief Конструктор для построения плоской матрицы упругости изотропного тела
     * @param E Модуль Юнга
     * @param nu Коэффициент Пуассона
     * @param isPlaneStress если true, то матрица строится для плоского напряженного состояния, если false - для плоского деформированного состояния
     */
    ElasticMatrix(double E, double nu, bool isPlaneStress);
    /**
     * @brief Конструктор копирования
     * @param ematrix Экземпляр объекта для копирования
     */
    ElasticMatrix(const ElasticMatrix &ematrix);
    /**
     * @brief Метод возвращает матрицу упругости
     * @return Матрица упругости
     */
    DoubleMatrix D() const;
    /**
     * @brief Установить значение матрицы упругости
     * @param D Матрица упругости
     */
    void setD(const DoubleMatrix &D);

private:
    DoubleMatrix D_; //!< Матрица упругости
};

#endif // ELASTICMATRIX_H
