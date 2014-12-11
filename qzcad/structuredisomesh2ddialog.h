/**
  * @author Сергей Чопоров
  * @date 07/12/2014
  * @version 1.0.1
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef STRUCTUREDISOMESH2DDIALOG_H
#define STRUCTUREDISOMESH2DDIALOG_H

#include <QDialog>

namespace Ui {
class StructuredIsoMesh2DDialog;
}
/**
 * @brief Класс формы для введения параметров изопараметрической сетки
 */
class StructuredIsoMesh2DDialog : public QDialog
{
    Q_OBJECT
    
public:
    /**
     * @brief Коснтруктор создает форму с параметрами полей по умолчанию
     * @param parent Указатель на родительский виджет
     */
    explicit StructuredIsoMesh2DDialog(QWidget *parent = 0);
    /**
     * @brief Конструктор создает форму с инициализацией полей заданными значениями
     * @param parent Указатель на родительский виджет
     * @param x0 Абсцисса точки 0
     * @param y0 Ордината точки 0
     * @param x1 Абсцисса точки 1
     * @param y1 Ордината точки 1
     * @param x2 Абсцисса точки 2
     * @param y2 Ордината точки 2
     * @param x3 Абсцисса точки 3
     * @param y3 Ордината точки 3
     * @param xiCount Колчество узлов сетки вдоль первого напраления
     * @param etaCount Количество узлов сетки вдоль второго направления
     */
    explicit StructuredIsoMesh2DDialog(QWidget *parent, const double &x0, const double &y0, const double &x1, const double &y1, const double &x2, const double &y2, const double &x3, const double &y3, const int &xiCount, const int &etaCount);
    ~StructuredIsoMesh2DDialog();
    /**
     * @brief Абсцисса точки 0
     * @return Значение абсциссы точки 0
     */
    double x0() const;
    /**
     * @brief Абсцисса точки 1
     * @return Значение абсциссы точки 1
     */
    double x1() const;
    /**
     * @brief Абсцисса точки 2
     * @return Значение абсциссы точки 2
     */
    double x2() const;
    /**
     * @brief Абсцисса точки 3
     * @return Значение абсциссы точки 3
     */
    double x3() const;
    /**
     * @brief Ордината точки 0
     * @return Значение ординаты точки 0
     */
    double y0() const;
    /**
     * @brief Ордината точки 1
     * @return Значение ординаты точки 1
     */
    double y1() const;
    /**
     * @brief Ордината точки 2
     * @return Значение ординаты точки 2
     */
    double y2() const;
    /**
     * @brief Ордината точки 3
     * @return Значение ординаты точки 3
     */
    double y3() const;
    /**
     * @brief x Значение ординаты i-й точки
     * @param i Номер точки (0 <= i < 4)
     * @return Значение ординаты точки с номером i
     */
    double x(int i) const;
    /**
     * @brief y Значение абсциссы i-й точки
     * @param i Номер точки (0 <= i < 4)
     * @return Значение абсциссы точки с номером i
     */
    double y(int i) const;
    /**
     * @brief xiCount Количество узлов вдоль первого направления
     * @return Количество узлов вдоль первого направления
     */ 
    int xiCount() const;
    /**
     * @brief etaCount Количество узлов вдоль второго направления
     * @return Количество узлов вдоль второго направления
     */
    int etaCount() const;
private:
    Ui::StructuredIsoMesh2DDialog *ui;
};

#endif // STRUCTUREDISOMESH2DDIALOG_H
