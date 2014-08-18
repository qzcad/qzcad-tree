#include "boundaryconditionswidget.h"
#include "ui_boundaryconditionswidget.h"
#include <QMessageBox>

BoundaryConditionsWidget::BoundaryConditionsWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::BoundaryConditionsWidget)
{
    ui->setupUi(this);
    ui->boundaryConditionsTable->setColumnWidth(1, 60);
    ui->boundaryConditionsTable->setColumnWidth(2, 60);
    ui->boundaryConditionsTable->setColumnWidth(3, 60);
}

BoundaryConditionsWidget::~BoundaryConditionsWidget()
{
    delete ui;
}

void BoundaryConditionsWidget::on_addCondition_clicked()
{
    ui->boundaryConditionsTable->insertRow(ui->boundaryConditionsTable->rowCount());
    QTableWidgetItem *uItem = new QTableWidgetItem("0.0");
    uItem->setCheckState(Qt::Checked);
    ui->boundaryConditionsTable->setItem(ui->boundaryConditionsTable->rowCount() - 1, 1, uItem);
    QTableWidgetItem *vItem = new QTableWidgetItem("0.0");
    vItem->setCheckState(Qt::Checked);
    ui->boundaryConditionsTable->setItem(ui->boundaryConditionsTable->rowCount() - 1, 2, vItem);
    QTableWidgetItem *wItem = new QTableWidgetItem("0.0");
    wItem->setCheckState(Qt::Checked);
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
