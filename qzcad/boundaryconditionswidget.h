#ifndef BOUNDARYCONDITIONSWIDGET_H
#define BOUNDARYCONDITIONSWIDGET_H

#include <QWidget>

namespace Ui {
class BoundaryConditionsWidget;
}
/**
 * @brief Виджет для ввода граничных условий и нагрузок
 */
class BoundaryConditionsWidget : public QWidget
{
    Q_OBJECT

public:
    explicit BoundaryConditionsWidget(QWidget *parent = 0);
    ~BoundaryConditionsWidget();
    /**
     * @brief Количество введенных условий
     * @return
     */
    int conditionsCount();
    /**
     * @brief Критерий отбора точек, попдающих под действие условия
     * @param i Номер условия
     * @return Текст с кодом функции критерия
     */
    QString condition(int i);
    /**
     * @brief Действие (в общем случае функция) в первом напралении
     * @param i Номер условия
     * @return Текст с кодом функции для действия в первом направлении
     */
    QString u(int i);
    /**
     * @brief Применять ли условие в первом направлении
     * @param i Номер условия
     * @return true, если условие применяется
     */
    bool isU(int i);
    /**
     * @brief Действие (в общем случае функция) во втором напралении
     * @param i Номер условия
     * @return Текст с кодом функции для действия во втором направлении
     */
    QString v(int i);
    /**
     * @brief Применять ли условие во втором направлении
     * @param i Номер условия
     * @return true, если условие применяется
     */
    bool isV(int i);
    /**
     * @brief Действие (в общем случае функция) в третьем напралении
     * @param i Номер условия
     * @return Текст с кодом функции для действия в третьем направлении
     */
    QString w(int i);
    /**
     * @brief Применять ли условие в третьем направлении
     * @param i Номер условия
     * @return true, если условие применяется
     */
    bool isW(int i);
    /**
     * @brief Установить режим работы: сила/закрепление
     * @param mode true -- ввод сил, false -- ввод граничных условий
     */
    void setIsForceMode(bool mode);
    /**
     * @brief Получить текущую размерность
     * @return true, если ввод для трехмерного пространства
     */
    bool getIs3d() const;
    /**
     * @brief Установить размерность пространтсва
     * @param value true, если ввод для трехмерного пространства
     */
    void setIs3d(bool value);
    /**
     * @brief Индекс типа нагрузки
     * @param i Номер нагрузки
     * @return Индекс типа нагрузки
     */
    int forceTypeIndex(int i);

private slots:
    void on_addCondition_clicked();

    void on_delCondition_clicked();

    void on_saveConditions_clicked();

    void on_loadConditions_clicked();

private:
    Ui::BoundaryConditionsWidget *ui;
    bool isForceMode;
    bool is3d;
};

#endif // BOUNDARYCONDITIONSWIDGET_H
