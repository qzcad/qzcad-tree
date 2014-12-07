/**
  * @author Сергей Чопоров
  * @date 07/12/2014
  * @version 1.0.1
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef POINTEDITORDIALOG_H
#define POINTEDITORDIALOG_H

#include <QDialog>

namespace Ui {
class PointEditorDialog;
}
/**
 * @brief Класс PointEditorDialog - редактор координат точки
 */
class PointEditorDialog : public QDialog
{
    Q_OBJECT
    
public:
    /**
     * @brief PointEditorDialog Конструктор редактора координат одномерной точки
     * @param parent Указатель на родительский виджет
     * @param x Значение ординаты
     */
    explicit PointEditorDialog(QWidget *parent, const double &x);
    /**
     * @brief PointEditorDialog Конструктор редактора двумерной точки
     * @param parent Указатель на родительский виджет
     * @param x Ордината
     * @param y Абсцисса
     */
    explicit PointEditorDialog(QWidget *parent, const double &x, const double &y);
    /**
     * @brief PointEditorDialog Конструктор редактора трехмерной точки
     * @param parent Указатель на родительский виджет
     * @param x Ордината
     * @param y Абсцисса
     * @param z Аппликата
     */
    explicit PointEditorDialog(QWidget *parent, const double &x, const double &y, const double &z);
    ~PointEditorDialog();
    /**
     * @brief xValue Ордината
     * @return Значение ординаты
     */
    double xValue();
    /**
     * @brief yValue Абсцисса
     * @return Значение абсциссы
     */
    double yValue();
    /**
     * @brief zValue Аппликата
     * @return Значение аппликаты
     */
    double zValue();
private:
    Ui::PointEditorDialog *ui;
};

#endif // POINTEDITORDIALOG_H
