#ifndef BOUNDARYCONDITIONSWIDGET_H
#define BOUNDARYCONDITIONSWIDGET_H

#include <QWidget>

namespace Ui {
class BoundaryConditionsWidget;
}

class BoundaryConditionsWidget : public QWidget
{
    Q_OBJECT

public:
    explicit BoundaryConditionsWidget(QWidget *parent = 0);
    ~BoundaryConditionsWidget();

private slots:
    void on_addCondition_clicked();

    void on_delCondition_clicked();

private:
    Ui::BoundaryConditionsWidget *ui;
};

#endif // BOUNDARYCONDITIONSWIDGET_H
