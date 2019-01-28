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
    ui->vectorScale->hide();
    ui->showInitialFrames->hide();
}

GLControlWidget::~GLControlWidget()
{
    delete ui;
}

GLMeshPicture *GLControlWidget::getGlMeshPicture()
{
    return ui->picture;
}

void GLControlWidget::activateDoubleBufferGL(bool activate)
{
    ui->picture->activateDoubleBufferGL(activate);
}

void GLControlWidget::activateTwoSideLightModel(bool activate)
{
    ui->picture->activateTwoSideLightModel(activate);
}

void GLControlWidget::activateSliceX(bool activate)
{
    ui->picture->activateSliceX(activate);
}

void GLControlWidget::activateSliceY(bool activate)
{
    ui->picture->activateSliceY(activate);
}

void GLControlWidget::activateSliceZ(bool activate)
{
    ui->picture->activateSliceZ(activate);
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
