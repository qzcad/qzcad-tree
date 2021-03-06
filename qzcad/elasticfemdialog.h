/**
  * @author Сергей Чопоров
  * @date 07/12/2014
  * @version 1.0.1
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef ELASTICFEMDIALOG_H
#define ELASTICFEMDIALOG_H

#include <QDialog>

namespace Ui {
class ElasticFemDialog;
}
/**
 * @brief Диалог для ввода параметров упругого расчета
 */
class ElasticFemDialog : public QDialog
{
    Q_OBJECT

public:
    explicit ElasticFemDialog(QWidget *parent = 0);
    ~ElasticFemDialog();

    void set2dMode(bool is2dMode = true);

    /**
     * @brief Изотропный материал
     * @return true, если материал изотропный
     */
    bool isIsotropy();
    /**
     * @brief Ортропный матриал (частный случай: E, nu, G)
     * @return true, если материал ортотропный
     */
    bool isOrthotropyENuG();
    /**
     * @brief Ортотропный материал (общий случай)
     * @return true, если материал ортотропный
     */
    bool isOrthotropy();
    /**
     * @brief Модуль Юнга (изотропный материал)
     * @return модуль Юнга для изотропного материала
     */
    double e();
    /**
     * @brief Модуль Юнга (ортотропный материал)
     * @return модуль Юнга для первого направления
     */
    double e1();
    /**
     * @brief Модуль Юнга (ортотропный материал)
     * @return модуль Юнга для второго направления
     */
    double e2();
    /**
     * @brief Модуль Юнга (ортотропный материал)
     * @return модуль Юнга для третьего направления
     */
    double e3();
    /**
     * @brief Коэффициент Пуассона (изотропный материал)
     * @return коэффициент Пуассона
     */
    double nu();
    /**
     * @brief Коэффициент Пуассона (ортотропный материал)
     * @return коэффициент Пуассона в плоскости 1-2
     */
    double nu12();
    /**
     * @brief Коэффициент Пуассона (ортотропный материал)
     * @return коэффициент Пуассона в плоскости 1-3
     */
    double nu13();
    /**
     * @brief Коэффициент Пуассона (ортотропный материал)
     * @return коэффициент Пуассона в плоскости 2-3
     */
    double nu23();
    /**
     * @brief Модуль сдвига (изотропный материал)
     * @return модуль сдвига
     */
    double g();
    /**
     * @brief Модуль сдвига (ортотропный материал)
     * @return Модуль сдвига в плоскости 1-2
     */
    double g12();
    /**
     * @brief Модуль сдвига (ортотропный материал)
     * @return Модуль сдвига в плоскости 1-3
     */
    double g13();
    /**
     * @brief Модуль сдвига (ортотропный материал)
     * @return Модуль сдвига в плоскости 2-3
     */
    double g23();
    /// Граничные условия
    /// @{
    /**
     * @brief Количество граничных условий
     * @return Количество граничных условий
     */
    int boundaryCount();
    /**
     * @brief Критерий отбора точек, попдающих под действие граничного условия
     * @param i Номер условия
     * @return Текст с кодом функции критерия
     */
    QString boundaryCondition(int i);
    /**
     * @brief Действие (в общем случае функция) в первом напралении
     * @param i Номер условия
     * @return Текст с кодом функции для действия в первом направлении
     */
    QString boundaryU(int i);
    /**
     * @brief Применять ли условие в первом направлении
     * @param i Номер условия
     * @return true, если условие применяется
     */
    bool boundaryIsU(int i);
    /**
     * @brief Действие (в общем случае функция) во втором напралении
     * @param i Номер условия
     * @return Текст с кодом функции для действия во втором направлении
     */
    QString boundaryV(int i);
    /**
     * @brief Применять ли условие во втором направлении
     * @param i Номер условия
     * @return true, если условие применяется
     */
    bool boundaryIsV(int i);
    /**
     * @brief Действие (в общем случае функция) в третьем напралении
     * @param i Номер условия
     * @return Текст с кодом функции для действия в третьем направлении
     */
    QString boundaryW(int i);
    /**
     * @brief Применять ли условие в третьем направлении
     * @param i Номер условия
     * @return true, если условие применяется
     */
    bool boundaryIsW(int i);
    /// @}
    /// Нагрузки
    /// @{
    /**
     * @brief Количество приложенных нагрузок
     * @return Количество приложенных нагрузок
     */
    int forcesCount();
    /**
     * @brief Критерий отбора узлов, попдающих под действие нагрузки
     * @param i Номер условия
     * @return Текст с кодом функции критерия
     */
    QString forceCondition(int i);
    /**
     * @brief Действие (в общем случае функция) в первом напралении
     * @param i Номер условия
     * @return Текст с кодом функции для действия в первом направлении
     */
    QString forceU(int i);
    /**
     * @brief Действие (в общем случае функция) во втором напралении
     * @param i Номер условия
     * @return Текст с кодом функции для действия во втором направлении
     */
    QString forceV(int i);
    /**
     * @brief Действие (в общем случае функция) в третьем напралении
     * @param i Номер условия
     * @return Текст с кодом функции для действия в третьем направлении
     */
    QString forceW(int i);
    /**
     * @brief Индекс типа нагрузки
     * @param i Номер нагрузки
     * @return Индекс типа нагрузки
     */
    int forceTypeIndex(int i);
    /// @}
private:
    Ui::ElasticFemDialog *ui;
};

#endif // ELASTICFEMDIALOG_H
