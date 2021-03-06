#include "exportmeshdialog.h"
#include "ui_exportmeshdialog.h"

ExportMeshDialog::ExportMeshDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ExportMeshDialog)
{
    ui->setupUi(this);
}

ExportMeshDialog::~ExportMeshDialog()
{
    delete ui;
}

bool ExportMeshDialog::isLayers()
{
    return ui->isLayers->isChecked();
}

