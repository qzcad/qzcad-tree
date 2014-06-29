/**
  * @author Сергей Чопоров
  * @date 19/06/2014
  * @version 1.0.1
  */
#ifndef NODE3D_H
#define NODE3D_H
#include "point3d.h"
#include "adjacentset.h"
#include "nodetype.h"

namespace msh
{
/**
  * @brief Узел в 3D
  */
typedef struct Node3D
{
    Point3D point; //!< Координаты узла
    AdjacentSet adjacent; //!< Список смежных элементов
    NodeType type; //!< Тип узла
} Node3D;
}
#endif // NODE3D_H
