/**
  * @author Сергей Чопоров
  * @date 25/02/2014
  * @version 1.0.0
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef NODETYPE_H
#define NODETYPE_H

namespace msh
{
/**
  * @brief Тип узла: внутренний, граничный, внешний или неопределенный (два последних необходимы для теоретической полноты)
  */
typedef enum { INNER = 1, BORDER = 0, OUTER = -1, CHARACTER = 2, UNDEFINED = -2 } NodeType;
}
#endif // NODETYPE_H
