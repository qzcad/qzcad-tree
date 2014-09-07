#include "elasticconstatntswidget.h"
#include "ui_elasticconstatntswidget.h"

ElasticConstatntsWidget::ElasticConstatntsWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ElasticConstatntsWidget)
{
    ui->setupUi(this);
    ui->e1->hide();
    ui->e1Label->hide();
    ui->e2->hide();
    ui->e2Label->hide();
    ui->e3->hide();
    ui->e3Label->hide();
    ui->nu12->hide();
    ui->nu12Label->hide();
    ui->nu13->hide();
    ui->nu13Label->hide();
    ui->nu23->hide();
    ui->nu23Label->hide();
    ui->shearGroupBox->hide();
    ui->g12->hide();
    ui->g12Label->hide();
    ui->g13->hide();
    ui->g13Label->hide();
    ui->g23->hide();
    ui->g23Label->hide();
    is3d = true;
}

ElasticConstatntsWidget::~ElasticConstatntsWidget()
{
    delete ui;
}
bool ElasticConstatntsWidget::getIs3d() const
{
    return is3d;
}

void ElasticConstatntsWidget::setIs3d(bool value)
{
    is3d = value;
    if (!is3d)
    {
        ui->e3->setEnabled(false);
        ui->e3Label->setEnabled(false);
        ui->nu13->setEnabled(false);
        ui->nu13Label->setEnabled(false);
        ui->nu23->setEnabled(false);
        ui->nu23Label->setEnabled(false);
        ui->g13->setEnabled(false);
        ui->g13Label->setEnabled(false);
        ui->g23->setEnabled(false);
        ui->g23Label->setEnabled(false);
    }
}

bool ElasticConstatntsWidget::isIsotropy()
{
    return ui->isIsotropy->isChecked();
}

bool ElasticConstatntsWidget::isOrthotropyENuG()
{
    return ui->isOrthotropyENuG->isChecked();
}

bool ElasticConstatntsWidget::isOrthotropy()
{
    return ui->isOrthotropy->isChecked();
}

double ElasticConstatntsWidget::e()
{
    return ui->e->value();
}

double ElasticConstatntsWidget::e1()
{
    return ui->e1->value();
}

double ElasticConstatntsWidget::e2()
{
    return ui->e2->value();
}

double ElasticConstatntsWidget::e3()
{
    return ui->e3->value();
}

double ElasticConstatntsWidget::nu()
{
    return ui->nu->value();
}

double ElasticConstatntsWidget::nu12()
{
    return ui->nu12->value();
}

double ElasticConstatntsWidget::nu13()
{
    return ui->nu13->value();
}

double ElasticConstatntsWidget::nu23()
{
    return ui->nu23->value();
}

double ElasticConstatntsWidget::g()
{
    return ui->g->value();
}

double ElasticConstatntsWidget::g12()
{
    return ui->g12->value();
}

double ElasticConstatntsWidget::g13()
{
    return ui->g13->value();
}

double ElasticConstatntsWidget::g23()
{
    return ui->g23->value();
}

