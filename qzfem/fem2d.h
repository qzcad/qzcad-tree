#ifndef FEM2D_H
#define FEM2D_H

#include <functional>

#include "fem.h"

#include "point2d.h"
using namespace msh;

#include "doublevector.h"
using namespace mtx;
/**
 * @brief Базовый класс для решения двумерных задач МКЭ
 */
class Fem2D : public Fem
{
public:
    /**
     * @brief Конструктор
     * @param mesh Указатель на сетку конечных элементов
     */
    Fem2D(Mesh *mesh, UInteger freedom, const std::list<FemCondition *> &conditions);
    /**
     * @brief Метод вычисляет матрицу плоского деформирвоанного состояния
     * @param E Модуль Юнга
     * @param nu Коэффициент Пуассона
     * @return Матрица 3x3 - матрица плоского деформированного состояния
     */
    static DoubleMatrix evalPlaneStrainMatrix(const double &E, const double &nu);
    /**
     * @brief Метод вычисляет матрицу плоского деформирвоанного состояния (случай независимого модуля сдвига)
     * @param E Модуль Юнга
     * @param nu Коэффициент Пуассона
     * @param G Модуль сдвига
     * @return Матрица 3x3 - матрица плоского деформированного состояния
     */
    static DoubleMatrix evalPlaneStrainMatrix(const double &E, const double &nu, const double &G);
    /**
     * @brief Метод вычисляет матрицу плоского напряженного состояния
     * @param E Модуль Юнга
     * @param nu Коэффициент Пуассона
     * @return Матрица 3x3 - матрица плоского напряженного состояния
     */
    static DoubleMatrix evalPlaneStressMatrix(const double &E, const double &nu);
    /**
     * @brief Метод вычисляет матрицу плоского напряженного состояния (случай независимого модуля сдвига)
     * @param E Модуль Юнга
     * @param nu Коэффициент Пуассона
     * @param G Модуль сдвига
     * @return Матрица 3x3 - матрица плоского напряженного состояния
     */
    static DoubleMatrix evalPlaneStressMatrix(const double &E, const double &nu, const double &G);
    /**
     * @brief Метод вычисляет матрицу плоского напряженного состояния ламината
     * @param E1 Модуль Юнга xx
     * @param E2 Модуль Юнга yy
     * @param nu12 Коэффициент Пуассона xy
     * @param nu21 Коэффициент Пуассона yx
     * @param G12 Модуль сдвига xy
     * @param theta Угол намотки
     * @return Матрица 3x3 - матрица плоского напряженного состояния ламината
     */
    static DoubleMatrix evalLaminaStressMatrix(const double &E1, const double &E2, const double &nu12, const double &G12, const double &theta);
    /**
     * @brief Метод вычисляет матрицу сдвига ламината
     * @param G13 Модуль свига xz
     * @param G23 Модуль сдвига yz
     * @param theta Угол намотки
     * @return Матрица 2x2 - матрица сдвига ламината
     */
    static DoubleMatrix evalLaminaShearMatrix(const double &G13, const double &G23, const double &theta);

protected:
    /**
     * @brief Метод для построения значений функций формы билинейного четырехугольного элемента
     * @param xi Значение параметра первого направления местной системы координат
     * @param eta Значения параметра второго направления местной системы координат
     * @param x Массив x-координат узлов
     * @param y Массив y-координат узлов
     * @param N Значения функций формы
     * @param dNdX Значения x-производной функций формы
     * @param dNdY Значения y-производной функций формы
     * @return Якобиан преобразования в местную систему координат
     */
    double isoQuad4(const double &xi, const double &eta, const DoubleVector &x, const DoubleVector &y,
                    DoubleVector &N, DoubleVector &dNdX, DoubleVector &dNdY);
    /**
     * @brief Метод для построения значений функций формы линейного треугольного элемента
     * @param xi Значение параметра первого направления местной системы координат
     * @param eta Значения параметра второго направления местной системы координат
     * @param x Массив x-координат узлов
     * @param y Массив y-координат узлов
     * @param N Значения функций формы
     * @param dNdX Значения x-производной функций формы
     * @param dNdY Значения y-производной функций формы
     * @return Якобиан преобразования в местную систему координат
     */
    double isoTriangle3(const double &xi, const double &eta, const DoubleVector &x, const DoubleVector &y,
                        DoubleVector &N, DoubleVector &dNdX, DoubleVector &dNdY);
};

#endif // FEM2D_H
