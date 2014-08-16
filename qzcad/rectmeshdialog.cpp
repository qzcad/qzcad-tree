#include "rectmeshdialog.h"
#include "ui_rectmeshdialog.h"

RectMeshDialog::RectMeshDialog(QWidget *parent, int dimension) :
    QDialog(parent),
    ui(new Ui::RectMeshDialog)
{
    ui->setupUi(this);

    if (dimension == 2)
    {
        ui->zCount->hide();
        ui->zCountLabel->hide();
        ui->zMin->hide();
        ui->zMinLabel->hide();
        ui->depth->hide();
        ui->depthLable->hide();
    }
}

RectMeshDialog::~RectMeshDialog()
{
    delete ui;
}

double RectMeshDialog::xMin()
{
    return ui->xMin->value();
}

double RectMeshDialog::yMin()
{
    return ui->yMin->value();
}

double RectMeshDialog::zMin()
{
    return ui->zMin->value();
}

double RectMeshDialog::rectWidth()
{
    return ui->width->value();
}

double RectMeshDialog::rectHeight()
{
    return ui->height->value();
}

double RectMeshDialog::rectDepth()
{
    return ui->depth->value();
}

int RectMeshDialog::xCount()
{
    return ui->xCount->value();
}

int RectMeshDialog::yCount()
{
    return ui->yCount->value();
}

int RectMeshDialog::zCount()
{
    return ui->zCount->value();
}
