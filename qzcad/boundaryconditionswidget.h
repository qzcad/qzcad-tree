#ifndef BOUNDARYCONDITIONSWIDGET_H
#define BOUNDARYCONDITIONSWIDGET_H

#include <QWidget>
#include "femcondition3d.h"

namespace Ui {
class BoundaryConditionsWidget;
}

class BoundaryConditionsWidget : public QWidget
{
    Q_OBJECT

public:
    explicit BoundaryConditionsWidget(QWidget *parent = 0);
    ~BoundaryConditionsWidget();
    int conditionsCount();
    FEMCondition3DPointer conditionPointer(int i);
    void setIsForceMode(bool mode);
private slots:
    void on_addCondition_clicked();

    void on_delCondition_clicked();

private:
    Ui::BoundaryConditionsWidget *ui;
    bool isForceMode;
};

#endif // BOUNDARYCONDITIONSWIDGET_H
