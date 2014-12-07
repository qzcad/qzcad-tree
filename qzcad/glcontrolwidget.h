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
    /**
     * @brief Получить текущий (выбранный) вектор значений
     * @return Вектор значений, который выбрал пользователь для визуализации
     * Если нет векторов в соответствующих массивах, будет возвращен пустой вектор.
     */
    NamedFloatingVector currentValuesVector();
public slots:
    void activateDoubleBufferGL(bool activate);
private slots:
    void on_mouseMode_currentIndexChanged(int index);

private:
    Ui::GLControlWidget *ui;
};

#endif // GLCONTROLWIDGET_H
