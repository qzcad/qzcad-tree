#include "glcontrolwidget.h"
#include "ui_glcontrolwidget.h"

GLControlWidget::GLControlWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::GLControlWidget)
{
    ui->setupUi(this);

    ui->labelTranslateStep->hide();
    ui->translateStep->hide();

    ui->labelXTranslation->hide();
    ui->labelYTranslation->hide();
    ui->labelZTranslation->hide();
    ui->xTranslation->hide();
    ui->yTranslation->hide();
    ui->zTranslation->hide();
}

GLControlWidget::~GLControlWidget()
{
    delete ui;
}

msh::MeshPointer GLControlWidget::getMesh()
{
    return ui->picture->getMesh();
}

void GLControlWidget::setMesh(msh::MeshPointer mesh)
{
    ui->picture->setMesh(mesh);
}

void GLControlWidget::resetMesh()
{
    ui->picture->resetMesh();
}

msh::MeshPointer GLControlWidget::releaseMesh()
{
    return ui->picture->releaseMesh();
}

void GLControlWidget::pushNodeValuesVector(const NamedFloatingVector &vector)
{
    ui->picture->pushNodeValuesVector(vector);
}

void GLControlWidget::clearNodeValues()
{
    ui->picture->clearNodeValues();
}

void GLControlWidget::pushElementValuesVector(const NamedFloatingVector &vector)
{
    ui->picture->pushElementValuesVector(vector);
}

void GLControlWidget::clearElementValues()
{
    ui->picture->clearElementValues();
}

void GLControlWidget::activateDoubleBufferGL(bool activate)
{
    ui->picture->activateDoubleBufferGL(activate);
}

void GLControlWidget::on_mouseMode_currentIndexChanged(int index)
{
    if (index == 1)
    {
        ui->labelTranslateStep->show();
        ui->translateStep->show();
        ui->labelXTranslation->show();
        ui->labelYTranslation->show();
        ui->labelZTranslation->show();
        ui->xTranslation->show();
        ui->yTranslation->show();
        ui->zTranslation->show();
    }
    else
    {
        ui->labelTranslateStep->hide();
        ui->translateStep->hide();
        ui->labelXTranslation->hide();
        ui->labelYTranslation->hide();
        ui->labelZTranslation->hide();
        ui->xTranslation->hide();
        ui->yTranslation->hide();
        ui->zTranslation->hide();
    }
}
