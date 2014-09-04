#include "boundaryconditionswidget.h"
#include "ui_boundaryconditionswidget.h"
#include <QMessageBox>
#include "qtscriptfemcondition3d.h"
#include "qtscriptforcecondition3d.h"

BoundaryConditionsWidget::BoundaryConditionsWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::BoundaryConditionsWidget)
{
    ui->setupUi(this);
    ui->boundaryConditionsTable->setColumnWidth(1, 100);
    ui->boundaryConditionsTable->setColumnWidth(2, 100);
    ui->boundaryConditionsTable->setColumnWidth(3, 100);
    isForceMode = false;
}

BoundaryConditionsWidget::~BoundaryConditionsWidget()
{
    delete ui;
}

int BoundaryConditionsWidget::conditionsCount()
{
    return ui->boundaryConditionsTable->rowCount();
}

FEMCondition3DPointer BoundaryConditionsWidget::conditionPointer(int i)
{
    FEMCondition3DPointer condPointer = NULL;
    if (isForceMode)
    {
        QString condition = ui->boundaryConditionsTable->item(i, 0)->text();
        QString u = ui->boundaryConditionsTable->item(i, 1)->text();
        QString v = ui->boundaryConditionsTable->item(i, 2)->text();
        QString w = ui->boundaryConditionsTable->item(i, 3)->text();
        condPointer = new QtScriptForceCondition3D(condition, u, v, w);
    }
    else
    {
        QString func = ui->boundaryConditionsTable->item(i, 0)->text();
        bool isU = ui->boundaryConditionsTable->item(i, 1)->checkState() == Qt::CheckState::Checked;
        double u = ui->boundaryConditionsTable->item(i, 1)->text().toDouble();
        bool isV = ui->boundaryConditionsTable->item(i, 2)->checkState() == Qt::CheckState::Checked;
        double v = ui->boundaryConditionsTable->item(i, 2)->text().toDouble();
        bool isW = ui->boundaryConditionsTable->item(i, 3)->checkState() == Qt::CheckState::Checked;
        double w = ui->boundaryConditionsTable->item(i, 3)->text().toDouble();
        condPointer = new QtScriptFemCondition3D(func, isU, u, isV, v, isW, w);
    }

    return condPointer;
}

void BoundaryConditionsWidget::setIsForceMode(bool mode)
{
    isForceMode = mode;
}

void BoundaryConditionsWidget::on_addCondition_clicked()
{
    ui->boundaryConditionsTable->insertRow(ui->boundaryConditionsTable->rowCount());
    QTableWidgetItem *uItem = new QTableWidgetItem("0.0");
    if (!isForceMode) uItem->setCheckState(Qt::Checked);
    ui->boundaryConditionsTable->setItem(ui->boundaryConditionsTable->rowCount() - 1, 1, uItem);
    QTableWidgetItem *vItem = new QTableWidgetItem("0.0");
    if (!isForceMode) vItem->setCheckState(Qt::Checked);
    ui->boundaryConditionsTable->setItem(ui->boundaryConditionsTable->rowCount() - 1, 2, vItem);
    QTableWidgetItem *wItem = new QTableWidgetItem("0.0");
    if (!isForceMode) wItem->setCheckState(Qt::Checked);
    ui->boundaryConditionsTable->setItem(ui->boundaryConditionsTable->rowCount() - 1, 3, wItem);

}

void BoundaryConditionsWidget::on_delCondition_clicked()
{
    QString question = tr("Удалить выбранные строки?");
    if (QMessageBox::question(this, question, question, QMessageBox::Yes | QMessageBox::No) == QMessageBox::Yes)
    {
        QModelIndexList indexes = ui->boundaryConditionsTable->selectionModel()->selection().indexes();
        for (int i = 0; i < indexes.count(); ++i)
        {
            QModelIndex index = indexes.at(i);
            ui->boundaryConditionsTable->removeRow(index.row());
        }
    }
}
