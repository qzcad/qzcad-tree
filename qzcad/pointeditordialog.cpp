#include "pointeditordialog.h"
#include "ui_pointeditordialog.h"

PointEditorDialog::PointEditorDialog(QWidget *parent, const double &x) :
    QDialog(parent),
    ui(new Ui::PointEditorDialog)
{
    ui->setupUi(this);

    ui->x->setValue(x);

    ui->y->hide();
    ui->yLabel->hide();
    ui->z->hide();
    ui->zLabel->hide();
}

PointEditorDialog::PointEditorDialog(QWidget *parent, const double &x, const double &y) :
    QDialog(parent),
    ui(new Ui::PointEditorDialog)
{
    ui->setupUi(this);

    ui->x->setValue(x);
    ui->y->setValue(y);

    ui->z->hide();
    ui->zLabel->hide();
}

PointEditorDialog::PointEditorDialog(QWidget *parent, const double &x, const double &y, const double &z) :
    QDialog(parent),
    ui(new Ui::PointEditorDialog)
{
    ui->setupUi(this);

    ui->x->setValue(x);
    ui->y->setValue(y);
    ui->z->setValue(z);
}

PointEditorDialog::~PointEditorDialog()
{
    delete ui;
}

double PointEditorDialog::xValue()
{
    return ui->x->value();
}

double PointEditorDialog::yValue()
{
    return ui->y->value();
}

double PointEditorDialog::zValue()
{
    return ui->z->value();
}
