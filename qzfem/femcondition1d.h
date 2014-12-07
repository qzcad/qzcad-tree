/**
  * @author Сергей Чопоров
  * @date 14/08/2014
  * @version 1.0.0
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  * */
#ifndef FEMCONDITION1D_H
#define FEMCONDITION1D_H
#include "femcondition.h"
/**
 * @brief Граничные условия для одномерного случая (виртуальный класс)
 */
class FEMCondition1D : public FEMCondition
{
public:
    /**
     * @brief Виртуальный деструктор
     */
    virtual ~FEMCondition1D() {}
    /**
     * @brief Значение граничного условия в первом направлении
     * @return Значение граничного условия в первом направлении
     */
    virtual double u() = 0;
    /**
     * @brief Проверить действие условия в первом направлении
     * @return true, если действует
     */
    virtual bool isU() = 0;
};
/**
 * @brief Указатель на одномерные граничные условия
 */
typedef FEMCondition1D* FEMCondition1DPointer;
#endif // FEMCONDITION1D_H
