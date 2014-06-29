#include "baryquadsdialog.h"
#include "ui_baryquadsdialog.h"

BaryQuadsDialog::BaryQuadsDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::BaryQuadsDialog)
{
    ui->setupUi(this);
}

BaryQuadsDialog::BaryQuadsDialog(QWidget *parent, const double &x0, const double &y0, const double &x1, const double &y1, const double &x2, const double &y2, const int &nodesCount) :
    QDialog(parent),
    ui(new Ui::BaryQuadsDialog)
{
    ui->setupUi(this);
    ui->x0->setValue(x0);
    ui->y0->setValue(y0);
    ui->x1->setValue(x1);
    ui->y1->setValue(y1);
    ui->x2->setValue(x2);
    ui->y2->setValue(y2);
    ui->nodesCount->setValue(nodesCount);
}

BaryQuadsDialog::~BaryQuadsDialog()
{
    delete ui;
}

double BaryQuadsDialog::x0() const
{
    return ui->x0->value();
}

double BaryQuadsDialog::y0() const
{
    return ui->y0->value();
}

double BaryQuadsDialog::x1() const
{
    return ui->x1->value();
}

double BaryQuadsDialog::y1() const
{
    return ui->y1->value();
}

double BaryQuadsDialog::x2() const
{
    return ui->x2->value();
}

double BaryQuadsDialog::y2() const
{
    return ui->y2->value();
}

double BaryQuadsDialog::x(const int &i) const
{
    switch (i)
    {
    case 0:
        return x0();
    case 1:
        return x1();
    case 2:
        return x2();
    default:
        return x0();
    }
}

double BaryQuadsDialog::y(const int &i) const
{
    switch (i)
    {
    case 0:
        return y0();
    case 1:
        return y1();
    case 2:
        return y2();
    default:
        return y0();
    }
}

int BaryQuadsDialog::nodesCount() const
{
    return ui->nodesCount->value();
}
