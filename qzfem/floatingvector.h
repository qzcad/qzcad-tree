/**
  * @author Сергей Чопоров
  * @date 01/08/2014
  * @version 1.0.0
  */
#ifndef FLOATINGVECTOR_H
#define FLOATINGVECTOR_H

#include <boost/numeric/ublas/vector.hpp>

namespace ublas =  boost::numeric::ublas;
using namespace msh;
/**
 * @brief Вектор чисел с плавающей точкой (boost)
 */
typedef ublas::vector<double> FloatingVector;

#endif // FLOATINGVECTOR_H
