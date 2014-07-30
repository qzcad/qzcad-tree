/**
  * @author Сергей Чопоров
  * @date 30/07/2014
  * @version 1.0.1
  */
#ifndef GLOBALMATRIX_H
#define GLOBALMATRIX_H

#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector_of_vector.hpp>
#include "floating.h"

namespace ublas =  boost::numeric::ublas;
using namespace msh;
/**
 * @brief Глобальная матрица - вектор сжатых векторов
 * @see boos::numeric::ublas
 */
typedef ublas::generalized_vector_of_vector< Floating, ublas::row_major, ublas::vector<ublas::compressed_vector<Floating> > > GlobalMatrix;

#endif // GLOBALMATRIX_H
