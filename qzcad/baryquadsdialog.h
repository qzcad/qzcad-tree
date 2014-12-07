/**
  * @author Сергей Чопоров
  * @date 07/12/2014
  * @version 1.0.1
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef BARYQUADSDIALOG_H
#define BARYQUADSDIALOG_H

#include <QDialog>

namespace Ui {
class BaryQuadsDialog;
}
/**
 * @brief Диалог ввода параметров треугольной области, для которой необходимо построить дискретизацию
 */
class BaryQuadsDialog : public QDialog
{
    Q_OBJECT
    
public:
    /**
     * @brief BaryQuadsDialog Конструктор по умолчанию
     * @param parent Указатель на родительский виджет
     */
    explicit BaryQuadsDialog(QWidget *parent = 0);
    /**
     * @brief BaryQuadsDialog Конструктор для инициализации полей заданными значениями
     * @param parent Указатель на родительский виджет
     * @param x0 Ордината точки 0
     * @param y0 Абсцисса точки 0
     * @param x1 Ордината точки 1
     * @param y1 Абсцисса точки 1
     * @param x2 Ордината точки 2
     * @param y2 Адсцисса точки 2
     * @param nodesCount Количество узлов сетки на одну сторону треугольника
     */
    explicit BaryQuadsDialog(QWidget *parent, const double &x0, const double &y0, const double &x1, const double &y1, const double &x2, const double &y2, const int &nodesCount);
    ~BaryQuadsDialog();
    double x0() const;
    double y0() const;
    double x1() const;
    double y1() const;
    double x2() const;
    double y2() const;
    double x(const int &i) const;
    double y(const int &i) const;
    int nodesCount() const;
private:
    Ui::BaryQuadsDialog *ui;
};

#endif // BARYQUADSDIALOG_H
