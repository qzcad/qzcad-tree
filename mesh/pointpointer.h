/**
  * @author Сергей Чопоров
  * @version 1.0.0
  * @date 11/01/2014
  */
#ifndef POINTPOINTER_H
#define POINTPOINTER_H
#include <memory>
#include <point.h>

namespace msh
{
/**
 * @brief "Умный" указатель на точку (c++11)
 */
typedef std::shared_ptr<Point> PointPointer;
}

#endif // POINTPOINTER_H
