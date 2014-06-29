/**
  * @author Сергей Чопоров
  * @date 25/02/2014
  * @version 1.0.0
  */
#ifndef NODETYPE_H
#define NODETYPE_H

namespace msh
{
/**
  * @brief Тип узла: внутренний, граничный или внешний (не обходим для теоретической полноты)
  */
typedef enum { INNER = 1, BORDER = 0, OUTER = -1, CHARACTER = 2 } NodeType;
}
#endif // NODETYPE_H
