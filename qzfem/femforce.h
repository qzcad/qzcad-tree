/**
  * @author Сергей Чопоров
  * @date 28/11/2014
  * @version 1.0.0
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  * */
#ifndef FEMFORCE_H
#define FEMFORCE_H
/**
 * @brief Тип нагрузки
 */
typedef enum
{
    NODAL_FORCE = 0, //!< Узловая нагрузка
    SURFACE_FORCE, //!< Поверхностная нагрузка (давление)
    VOLUME_FORCE //!< Объемная сила
} ForceType;

/**
 * @brief Базовый абстрактный класс для передачи цсловия нагружения
 */
class FEMForce
{
public:
    /**
     * @brief Конструктор по умолчанию
     */
    FEMForce();
    /**
     * @brief Конструктор с установлением типа нагрузки
     * @param forceType Тип нагрузки
     * @see ForceType
     */
    FEMForce(ForceType forceType);
    /**
     * @brief Конструктор копирования
     * @param force Экземпляр объекта для создания копии
     */
    FEMForce(const FEMForce &force);
    /**
     * @brief Виртуальный деструктор
     */
    virtual ~FEMForce() {}
    /**
     * @brief Получить тип нагрузки
     * @return Тип нагрузки
     * @see ForceType
     */
    ForceType forceType() const;
    /**
     * @brief Установть тип нагрузки
     * @param forceType Новое значение типа нагрузки
     */
    void setForceType(const ForceType &forceType);

private:
    ForceType forceType_; //!< Тип нагрузки
};

#endif // FEMFORCE_H
