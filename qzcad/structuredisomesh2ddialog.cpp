#include "structuredisomesh2ddialog.h"
#include "ui_structuredisomesh2ddialog.h"

StructuredIsoMesh2DDialog::StructuredIsoMesh2DDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::StructuredIsoMesh2DDialog)
{
    ui->setupUi(this);
}

StructuredIsoMesh2DDialog::StructuredIsoMesh2DDialog(QWidget *parent, const double &x0, const double &y0, const double &x1, const double &y1, const double &x2, const double &y2, const double &x3, const double &y3, const int &xiCount, const int &etaCount) :
    QDialog(parent),
    ui(new Ui::StructuredIsoMesh2DDialog)
{
    ui->setupUi(this);
    ui->x0->setValue(x0);
    ui->y0->setValue(y0);
    ui->x1->setValue(x1);
    ui->y1->setValue(y1);
    ui->x2->setValue(x2);
    ui->y2->setValue(y2);
    ui->x3->setValue(x3);
    ui->y3->setValue(y3);
    ui->xiCount->setValue(xiCount);
    ui->etaCount->setValue(etaCount);
}

StructuredIsoMesh2DDialog::~StructuredIsoMesh2DDialog()
{
    delete ui;
}

double StructuredIsoMesh2DDialog::x0() const
{
    return ui->x0->value();
}

double StructuredIsoMesh2DDialog::x1() const
{
    return ui->x1->value();
}

double StructuredIsoMesh2DDialog::x2() const
{
    return ui->x2->value();
}

double StructuredIsoMesh2DDialog::x3() const
{
    return ui->x3->value();
}

double StructuredIsoMesh2DDialog::y0() const
{
    return ui->y0->value();
}

double StructuredIsoMesh2DDialog::y1() const
{
    return ui->y1->value();
}

double StructuredIsoMesh2DDialog::y2() const
{
    return ui->y2->value();
}

double StructuredIsoMesh2DDialog::y3() const
{
    return ui->y3->value();
}

double StructuredIsoMesh2DDialog::x(int i) const
{
    switch(i)
    {
    case 0:
        return ui->x0->value();
    case 1:
        return ui->x1->value();
    case 2:
        return ui->x2->value();
    case 3:
        return ui->x3->value();
    default:
        return ui->x0->value();
    }
}

double StructuredIsoMesh2DDialog::y(int i) const
{
    switch(i)
    {
    case 0:
        return ui->y0->value();
    case 1:
        return ui->y1->value();
    case 2:
        return ui->y2->value();
    case 3:
        return ui->y3->value();
    default:
        return ui->y0->value();
    }
}

int StructuredIsoMesh2DDialog::xiCount() const
{
    return ui->xiCount->value();
}

int StructuredIsoMesh2DDialog::etaCount() const
{
    return ui->etaCount->value();
}
