#ifndef GLCONTROLWIDGET_H
#define GLCONTROLWIDGET_H

#include <QWidget>
#include "meshpointer.h"

namespace Ui {
class GLControlWidget;
}

class GLControlWidget : public QWidget
{
    Q_OBJECT
    
public:
    explicit GLControlWidget(QWidget *parent = 0);
    ~GLControlWidget();
    msh::MeshPointer getMesh();
    void setMesh(msh::MeshPointer mesh);
    void resetMesh();
    msh::MeshPointer releaseMesh();
public slots:
    void activateDoubleBufferGL(bool activate);
private slots:
    void on_mouseMode_currentIndexChanged(int index);

private:
    Ui::GLControlWidget *ui;
};

#endif // GLCONTROLWIDGET_H
