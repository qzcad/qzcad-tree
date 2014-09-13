/**
  * @author Сергей Чопоров
  * @date 15/10/2012
  * @version 1.0.5
  */
#ifndef FLOATING_H
#define FLOATING_H

#include <float.h>

namespace msh
{
/**
 * @typedef Floating
 * Базовый тип для операций с плавающей точкой
 */
typedef double Floating;
/// Максимальное значение, которое может хнаить тип Floating
#undef FLOATING_MAX
#define FLOATING_MAX DBL_MAX
/// Минимальное значение, которое может хранить тип Floating
#undef FLOATING_MIN
#define FLOATING_MIN -DBL_MAX
}

#endif // FLOATING_H
