/**
  * @author Сергей Чопоров
  * @date 04/08/2014
  * @version 1.0.0
  */
#ifndef DETERMINANT_HPP
#define DETERMINANT_HPP

#include <boost/numeric/ublas/lu.hpp>
#include "floatingmatrix.h"

namespace ublas = boost::numeric::ublas;
/**
 * @brief Вычисление знака определителя
 * @param pm Матрица перестановок
 * @return Знак определителя
 */
int determinantSign(const ublas::permutation_matrix<std::size_t> &pm)
{
    int pm_sign = 1;
    std::size_t size = pm.size();
    for (std::size_t i = 0; i < size; ++i)
        if (i != pm(i))
            pm_sign *= -1; // swap_rows would swap a pair of rows here, so we change sign
    return pm_sign;
}
/**
 * @brief Функция для вычисления определителя матрицы
 * @param m Исходная матрица
 * @return Значение определителя
 */
double determinant(const FloatingMatrix &m)
{
    FloatingMatrix a(m); // рабочая копия матрицы
    ublas::permutation_matrix<std::size_t> pm(a.size1());
    double det = 1.0;
    if( ublas::lu_factorize(a, pm) )
    {
        det = 0.0;
    } else
    {
        for(FloatingMatrix::size_type i = 0; i < a.size1(); i++)
            det *= a(i, i); // multiply by elements on diagonal
        det = det * (double)determinantSign( pm );
    }
    return det;
}

#endif // DETERMINANT_HPP
