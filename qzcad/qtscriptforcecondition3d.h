#ifndef QTSCRIPTFORCECONDITION3D_H
#define QTSCRIPTFORCECONDITION3D_H

#include "femcondition3d.h"
#include <QtScript/QScriptEngine>

class QtScriptForceCondition3D: public FEMCondition3D
{
public:
    QtScriptForceCondition3D(QString condition, QString u, QString v, QString w);

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
    QScriptEngine engine_;
    QScriptValue object_;
    double x_;
    double y_;
    double z_;
};

#endif // QTSCRIPTFORCECONDITION3D_H
