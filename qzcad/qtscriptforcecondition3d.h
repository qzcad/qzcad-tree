/**
  * @author Сергей Чопоров
  * @date 24/08/2014
  * @version 1.0.1
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef QTSCRIPTFORCECONDITION3D_H
#define QTSCRIPTFORCECONDITION3D_H

#include "forcecondition3d.h"
#include "qzscriptengine.h"

class QtScriptForceCondition3D: public ForceCondition3D
{
public:
    QtScriptForceCondition3D(QString condition, QString u, QString v, QString w, ForceType fType);

    bool setCondition(QString condition, QString u, QString v, QString w);

    virtual bool isU();

    virtual double u();

    virtual bool isV();

    virtual double v();

    virtual bool isW();

    virtual double w();

    virtual bool isApplied(msh::PointPointer point);

    virtual ~QtScriptForceCondition3D(){}
private:
    QString condition_;
    QString u_;
    QString v_;
    QString w_;
    QZScriptEngine engine_;
    QScriptValue object_;
    double x_;
    double y_;
    double z_;
};

#endif // QTSCRIPTFORCECONDITION3D_H
