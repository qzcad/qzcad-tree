/**
  * @author Сергей Чопоров
  * @date 07/12/2014
  * @version 1.0.1
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef GLCONTROLWIDGET_H
#define GLCONTROLWIDGET_H

#include <QWidget>
#include "glmeshpicture.h"
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
    /**
     * @brief Получить указатель на экземпляр отрисовщика сеток
     * @return Укзатаель типа GLMeshPicture
     * @see GLMeshPicture
     */
    GLMeshPicture *getGlMeshPicture();
public slots:
    void activateDoubleBufferGL(bool activate);
    void activateTwoSideLightModel(bool activate);
    void activateSliceX(bool activate);
    void activateSliceY(bool activate);
    void activateSliceZ(bool activate);
private slots:
    void on_mouseMode_currentIndexChanged(int index);

private:
    Ui::GLControlWidget *ui;
};

#endif // GLCONTROLWIDGET_H
