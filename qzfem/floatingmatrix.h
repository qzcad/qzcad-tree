/**
  * @author Сергей Чопоров
  * @date 01/08/2014
  * @version 1.0.1
  */
#ifndef FLOATINGMATRIX_H
#define FLOATINGMATRIX_H

#include <boost/numeric/ublas/matrix.hpp>
#include "floating.h"

namespace ublas =  boost::numeric::ublas;
using namespace msh;
/**
 * @brief Матрица чисел с плавающей точкой (boost)
 */
typedef ublas::matrix<Floating> FloatingMatrix;

#endif // FLOATINGMATRIX_H
