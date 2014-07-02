#include "rotationbodymeshdialog.h"
#include "ui_rotationbodymeshdialog.h"

RotationBodyMeshDialog::RotationBodyMeshDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::RotationBodyMeshDialog)
{
    ui->setupUi(this);
    ui->labelRotationAngle->hide();
    ui->rotationAngle->hide();
}

RotationBodyMeshDialog::~RotationBodyMeshDialog()
{
    delete ui;
}

int RotationBodyMeshDialog::axe()
{
    return ui->axeComboBox->currentIndex();
}

double RotationBodyMeshDialog::radius()
{
    return ui->radius->value();
}

int RotationBodyMeshDialog::layersCount()
{
    return ui->layersCount->value();
}

bool RotationBodyMeshDialog::isClosedBody()
{
    return ui->isClosedBody->isChecked();
}

double RotationBodyMeshDialog::rotationAngle()
{
    return ui->rotationAngle->value();
}
