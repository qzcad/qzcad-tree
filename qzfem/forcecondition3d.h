/**
  * @author Сергей Чопоров
  * @date 28/11/2014
  * @version 1.0.0
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  * */
#ifndef FORCECONDITION3D_H
#define FORCECONDITION3D_H

#include "femcondition3d.h"
#include "femforce.h"

/**
 * @brief Виртуальный класс для задания нагрузок в трехмерном случае
 */
class ForceCondition3D : public FEMCondition3D, public FEMForce
{
public:
    /**
     * @brief Виртуальны деструктор
     */
    virtual ~ForceCondition3D() {}
};
/**
 * @brief Указатель на ля нагрузок в трехмерном случае
 */
typedef ForceCondition3D* ForceCondition3DPointer;
#endif // FORCECONDITION3D_H
