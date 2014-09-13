#ifndef GLCONTROLWIDGET_H
#define GLCONTROLWIDGET_H

#include <QWidget>
#include "meshpointer.h"
#include "namedfloatingvector.h"

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
    void pushNodeValuesVector(const NamedFloatingVector &vector);
    void clearNodeValues();
    /**
     * @brief Добавить вектор значений, определенных на элементе
     * @param vector Вектор значений определенных на элементе
     */
    void pushElementValuesVector(const NamedFloatingVector &vector);
    /**
     * @brief Очистить вектор значений, определенных на элементе
     */
    void clearElementValues();
public slots:
    void activateDoubleBufferGL(bool activate);
private slots:
    void on_mouseMode_currentIndexChanged(int index);

private:
    Ui::GLControlWidget *ui;
};

#endif // GLCONTROLWIDGET_H
