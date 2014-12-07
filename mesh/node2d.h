/**
  * @author Сергей Чопоров
  * @date 25/02/2014
  * @version 1.0.1
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef NODE2D_H
#define NODE2D_H
#include "point2d.h"
#include "adjacentset.h"
#include "nodetype.h"

namespace msh
{
/**
  * @brief Узел в 2D
  */
typedef struct Node2D
{
    Point2D point; //!< Координаты узла
    AdjacentSet adjacent; //!< Список смежных элементов
    NodeType type; //!< Тип узла
} Node2D;
}
#endif // NODE2D_H
