/**
  * @author Сергей Чопоров
  * @date 04/08/2014
  * @version 1.0.0
  */
#ifndef INVERTMATRIX_HPP
#define INVERTMATRIX_HPP
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>

namespace ublas = boost::numeric::ublas;

/**
 * @brief Функция для нахождения обратной матрицы (использует LU-разложение)
 * @param input Исходная матрица
 * @param inverse Обратная матрица
 */
template<class T>
bool invertMatrix (const ublas::matrix<T>& input, ublas::matrix<T>& inverse) {
    using namespace boost::numeric::ublas;
    typedef permutation_matrix<std::size_t> pmatrix;
    // рабочая копия исходной матрицы
    matrix<T> A(input);
    // матрица перестановок для LU-разложения
    pmatrix pm(A.size1());

    // LU-разложение
    int res = lu_factorize(A, pm);
    if( res != 0 ) return false;

    // создание единичной матрицы
    inverse.assign(ublas::identity_matrix<T>(A.size1()));

    // получение обратной матрицы
    lu_substitute(A, pm, inverse);

    return true;
}
#endif // INVERTMATRIX_HPP
