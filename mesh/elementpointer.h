/**
  * @author Сергей Чопоров
  * @version 1.0.0
  * @date 11/01/2014
  */
#ifndef ELEMENTPOINTER_H
#define ELEMENTPOINTER_H
#include <memory>
#include "element.h"
namespace msh
{
/**
 * @brief "Умный" указатель на элемент (c++11)
 */
typedef std::shared_ptr<Element> ElementPointer;
}
#endif // ELEMENTPOINTER_H
