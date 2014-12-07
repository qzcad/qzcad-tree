/**
  * @author Сергей Чопоров
  * @date 19/06/2014
  * @version 1.0.1
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef INTEGERVECTOR_H
#define INTEGERVECTOR_H

#include <vector>
#include "integer.h"

namespace msh {
/**
 * @brief Массив целых чисел (без знака)
 * @see UInteger
 */
typedef std::vector<UInteger> UIntegerVector;
}

#endif // INTEGERVECTOR_H
