#include "elasticfemdialog.h"
#include "ui_elasticfemdialog.h"

ElasticFemDialog::ElasticFemDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ElasticFemDialog)
{
    ui->setupUi(this);
}

ElasticFemDialog::~ElasticFemDialog()
{
    delete ui;
}
