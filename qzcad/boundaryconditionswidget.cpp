#include "boundaryconditionswidget.h"
#include "ui_boundaryconditionswidget.h"
#include <QMessageBox>
#include <QFile>
#include <QFileDialog>
#include <QTextStream>
#include <iostream>
#include <QDomDocument>

BoundaryConditionsWidget::BoundaryConditionsWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::BoundaryConditionsWidget)
{
    ui->setupUi(this);
    ui->boundaryConditionsTable->setColumnWidth(1, 100);
    ui->boundaryConditionsTable->setColumnWidth(2, 100);
    ui->boundaryConditionsTable->setColumnWidth(3, 100);
    isForceMode = false;
    is3d = true;
}

BoundaryConditionsWidget::~BoundaryConditionsWidget()
{
    delete ui;
}

int BoundaryConditionsWidget::conditionsCount()
{
    return ui->boundaryConditionsTable->rowCount();
}

QString BoundaryConditionsWidget::condition(int i)
{
    return ui->boundaryConditionsTable->item(i, 0)->text();
}

QString BoundaryConditionsWidget::u(int i)
{
    return ui->boundaryConditionsTable->item(i, 1)->text();
}

bool BoundaryConditionsWidget::isU(int i)
{
    return (ui->boundaryConditionsTable->item(i, 1)->checkState() == Qt::Checked);
}

QString BoundaryConditionsWidget::v(int i)
{
    return ui->boundaryConditionsTable->item(i, 2)->text();
}

bool BoundaryConditionsWidget::isV(int i)
{
    return (ui->boundaryConditionsTable->item(i, 2)->checkState() == Qt::Checked);
}

QString BoundaryConditionsWidget::w(int i)
{
    return ui->boundaryConditionsTable->item(i, 3)->text();
}

bool BoundaryConditionsWidget::isW(int i)
{
    return (ui->boundaryConditionsTable->item(i, 3)->checkState() == Qt::Checked);
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
bool BoundaryConditionsWidget::getIs3d() const
{
    return is3d;
}

void BoundaryConditionsWidget::setIs3d(bool value)
{
    is3d = value;
    if (!is3d)
    {
        ui->boundaryConditionsTable->setColumnHidden(3, true);
    }
}


void BoundaryConditionsWidget::on_saveConditions_clicked()
{
    QString fileName = QFileDialog::getSaveFileName(this, tr("Сохранить список условий"), "", tr("xml-файлы (*.xml)"));
    if(!fileName.isEmpty())
    {
        QFile file(fileName);
        if (!file.open(QFile::WriteOnly | QFile::Text))
        {
            std::cerr << "Ошибка: Невозможно открыть файл для записи " << qPrintable(fileName) << ": " << qPrintable(file.errorString()) << std::endl;
            return;
        }
        QDomDocument doc;
        QDomElement conditions = doc.createElement("conditions");
        conditions.setAttribute("isForceMode", isForceMode);
        for (int i = 0; i < ui->boundaryConditionsTable->rowCount(); i++)
        {
            QDomElement cond = doc.createElement("condition");
            cond.setAttribute("isApplied", condition(i));
            cond.setAttribute("u", u(i));
            cond.setAttribute("v", v(i));
            cond.setAttribute("w", w(i));
            if (!isForceMode)
            {
                cond.setAttribute("isU", isU(i));
                cond.setAttribute("isV", isV(i));
                cond.setAttribute("isW", isW(i));
            }
            conditions.appendChild(cond);
        }
        doc.appendChild(conditions);
        QTextStream out(&file);
        out << doc.toString();
        file.close();
    }
}

void BoundaryConditionsWidget::on_loadConditions_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Открыть список условий"), "", tr("xml-файлы (*.xml);;Любой файл (*)"));
    if(!fileName.isEmpty())
    {
        ui->boundaryConditionsTable->clear();
        while (ui->boundaryConditionsTable->rowCount() > 0)
        {
            ui->boundaryConditionsTable->removeRow(0);
        }

        QFile file(fileName);
        if (!file.open(QFile::ReadOnly | QFile::Text))
        {
            std::cerr << "Ошибка: Невозможно прочитать файл " << qPrintable(fileName) << ": " << qPrintable(file.errorString()) << std::endl;
            return;
        }
        QString errorStr;
        int errorLine;
        int errorColumn;
        QDomDocument doc;
        if (!doc.setContent(&file, false, &errorStr, &errorLine, &errorColumn))
        {
            std::cerr << "Ошибка: При разборе файла в сточке " << errorLine << ", столбце " << errorColumn << ": " << qPrintable(errorStr) << std::endl;
            return;
        }
        QDomElement root = doc.documentElement();
        if (root.tagName() != "conditions")
        {
            std::cerr << "Ошибка: Файл не является списком условий" << std::endl;
            return;
        }
        if (isForceMode != (root.attribute("isForceMode").toInt() == 1))
        {
            std::cerr << "Ошибка: ";
            if (!isForceMode)
            {
                std::cerr << "Файл не является списком граничных условий";
            }
            else
            {
                std::cerr << "Файл не является списком нагрузок";
            }
            std::cerr << std::endl;
            return;
        }
        QDomNode child = root.firstChild();
        while (!child.isNull())
        {
            QDomElement cond = child.toElement();
            if(cond.tagName() != "condition")
            {
                std::cerr << "Ошибка: Файл не является списком условий" << std::endl;
                return;
            }

            ui->boundaryConditionsTable->insertRow(ui->boundaryConditionsTable->rowCount());

            QTableWidgetItem *condItem = new QTableWidgetItem(cond.attribute("isApplied"));
            ui->boundaryConditionsTable->setItem(ui->boundaryConditionsTable->rowCount() - 1, 0, condItem);

            QTableWidgetItem *uItem = new QTableWidgetItem(cond.attribute("u"));
            if (!isForceMode)
            {
                bool isU = (cond.attribute("isU").toInt() == 1);
                uItem->setCheckState((isU == true) ? Qt::Checked : Qt::Unchecked);
            }
            ui->boundaryConditionsTable->setItem(ui->boundaryConditionsTable->rowCount() - 1, 1, uItem);

            QTableWidgetItem *vItem = new QTableWidgetItem(cond.attribute("v"));
            if (!isForceMode)
            {
                bool isV = (cond.attribute("isV").toInt() == 1);
                vItem->setCheckState((isV == true) ? Qt::Checked : Qt::Unchecked);
            }
            ui->boundaryConditionsTable->setItem(ui->boundaryConditionsTable->rowCount() - 1, 2, vItem);

            QTableWidgetItem *wItem = new QTableWidgetItem(cond.attribute("w"));
            if (!isForceMode)
            {
                bool isW = (cond.attribute("isW").toInt() == 1);
                wItem->setCheckState((isW == true) ? Qt::Checked : Qt::Unchecked);
            }
            ui->boundaryConditionsTable->setItem(ui->boundaryConditionsTable->rowCount() - 1, 3, wItem);

            child = child.nextSibling();
        }
    }
}
