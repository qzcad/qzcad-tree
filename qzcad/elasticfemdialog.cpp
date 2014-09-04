#include "elasticfemdialog.h"
#include "ui_elasticfemdialog.h"

ElasticFemDialog::ElasticFemDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ElasticFemDialog)
{
    ui->setupUi(this);
    ui->forces->setIsForceMode(true);
}

ElasticFemDialog::~ElasticFemDialog()
{
    delete ui;
}
