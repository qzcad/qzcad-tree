/**
  * @author Сергей Чопоров
  * @date 14/08/2014
  * @version 1.0.0
  * */
#ifndef FEMCONDITION2D_H
#define FEMCONDITION2D_H
#include "femcondition1d.h"
/**
 * @brief Граничные условия для двумерных задач (виртуальный класс)
 */
class FEMCondition2D : public FEMCondition1D
{
public:
    /**
     * @brief Виртуальный деструктор
     */
    virtual ~FEMCondition2D() {}
    /**
     * @brief Значение граничного условия во втором направлении
     * @return Значение граничного условия во втором направлении
     */
    virtual double v() = 0;
    /**
     * @brief Проверить действие условия во втором направлении
     * @return true, если действует во втором направлении
     */
    virtual bool isV() = 0;
};
/**
 * @brief Указатель на двумерные граничные условия
 */
typedef FEMCondition2D* FEMCondition2DPointer;
#endif // FEMCONDITION2D_H
