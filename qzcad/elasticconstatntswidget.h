/**
  * @author Сергей Чопоров
  * @date 04/09/2014
  * @version 1.0.0
  * */
#ifndef ELASTICCONSTATNTSWIDGET_H
#define ELASTICCONSTATNTSWIDGET_H

#include <QWidget>

namespace Ui {
class ElasticConstatntsWidget;
}
/**
 * @brief Класс виджета для ввода механических констант
 */
class ElasticConstatntsWidget : public QWidget
{
    Q_OBJECT

public:
    explicit ElasticConstatntsWidget(QWidget *parent = 0);
    ~ElasticConstatntsWidget();

    bool getIs3d() const;
    /**
     * @brief Изменить состояние флага трехмерного расчета
     * @param value новой значение флага
     */
    void setIs3d(bool value);
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

private:
    Ui::ElasticConstatntsWidget *ui;
    bool is3d;
};

#endif // ELASTICCONSTATNTSWIDGET_H
