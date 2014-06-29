/**
  * @author Сергей Чопоров
  * @version 1.0.0
  * @date 11/01/2014
  */
#ifndef MESHPOINTER_H
#define MESHPOINTER_H
#include <memory>
#include "mesh.h"

namespace msh
{
/**
 * @brief "Умный" указатель на сетку (c++11)
 */
typedef std::shared_ptr<Mesh> MeshPointer;
}
#endif // MESHPOINTER_H
