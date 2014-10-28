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

void ElasticFemDialog::set2dMode(bool is2dMode)
{
    if (is2dMode)
    {
        ui->constants->setIs3d(false);
        ui->boundaryConditions->setIs3d(false);
        ui->forces->setIs3d(false);
    }
    else
    {
        ui->constants->setIs3d(true);
        ui->boundaryConditions->setIs3d(true);
        ui->forces->setIs3d(true);
    }
}

bool ElasticFemDialog::isIsotropy()
{
    return ui->constants->isIsotropy();
}

bool ElasticFemDialog::isOrthotropyENuG()
{
    return ui->constants->isOrthotropyENuG();
}

bool ElasticFemDialog::isOrthotropy()
{
    return ui->constants->isOrthotropy();
}

double ElasticFemDialog::e()
{
    return ui->constants->e();
}

double ElasticFemDialog::e1()
{
    return ui->constants->e1();
}

double ElasticFemDialog::e2()
{
    return ui->constants->e2();
}

double ElasticFemDialog::e3()
{
    return ui->constants->e3();
}

double ElasticFemDialog::nu()
{
    return ui->constants->nu();
}

double ElasticFemDialog::nu12()
{
    return ui->constants->nu12();
}

double ElasticFemDialog::nu13()
{
    return ui->constants->nu13();
}

double ElasticFemDialog::nu23()
{
    return ui->constants->nu23();
}

double ElasticFemDialog::g()
{
    return ui->constants->g();
}

double ElasticFemDialog::g12()
{
    return ui->constants->g12();
}

double ElasticFemDialog::g13()
{
    return ui->constants->g13();
}

double ElasticFemDialog::g23()
{
    return ui->constants->g23();
}

int ElasticFemDialog::boundaryCount()
{
    return ui->boundaryConditions->conditionsCount();
}

QString ElasticFemDialog::boundaryCondition(int i)
{
    return ui->boundaryConditions->condition(i);
}

QString ElasticFemDialog::boundaryU(int i)
{
    return ui->boundaryConditions->u(i);
}

bool ElasticFemDialog::boundaryIsU(int i)
{
    return ui->boundaryConditions->isU(i);
}

QString ElasticFemDialog::boundaryV(int i)
{
    return ui->boundaryConditions->v(i);
}

bool ElasticFemDialog::boundaryIsV(int i)
{
    return ui->boundaryConditions->isV(i);
}

QString ElasticFemDialog::boundaryW(int i)
{
    return ui->boundaryConditions->w(i);
}

bool ElasticFemDialog::boundaryIsW(int i)
{
    return ui->boundaryConditions->isW(i);
}

int ElasticFemDialog::forcesCount()
{
    return ui->forces->conditionsCount();
}

QString ElasticFemDialog::forceCondition(int i)
{
    return ui->forces->condition(i);
}

QString ElasticFemDialog::forceU(int i)
{
    return ui->forces->u(i);
}

QString ElasticFemDialog::forceV(int i)
{
    return ui->forces->v(i);
}

QString ElasticFemDialog::forceW(int i)
{
    return ui->forces->w(i);
}

int ElasticFemDialog::forceTypeIndex(int i)
{
    return ui->forces->forceTypeIndex(i);
}
