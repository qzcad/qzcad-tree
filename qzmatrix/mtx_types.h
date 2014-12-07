/**
  * @author Сергей Чопоров
  * @date 23/09/2014
  * @version 1.0.0
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  * */
#ifndef MTX_TYPES_H
#define MTX_TYPES_H

namespace mtx {
/**
 * @brief Определение типа размера массива
 */
typedef unsigned long size_type;
/**
 * @brief Определение типа ссылок
 */
typedef double & reference;
/**
 * @brief Определение типа константных ссылок
 */
typedef const double & const_reference;
/**
 * @brief Определение типа указателей
 */
typedef double * pointer;
/**
 * @brief Определение типа константных указателей
 */
typedef const double * const_pointer;

}

#endif // MTX_TYPES_H
