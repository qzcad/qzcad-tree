#include <math.h>
#include <QObject>
#include <iostream>
#include "qpoint2d.h"
#include "qtscriptfunctions.h"

QScriptValue approx(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);
    double eps = 1.0E-6;
    if (context->argumentCount() != 2 && context->argumentCount() != 3)
    {
        return context->throwError(QObject::tr("approx() takes two or three arguments: approx(x, y) or approx(x, y, eps)"));
    }
    if (context->argumentCount() == 3)
    {
        if (!context->argument(2).isNumber())
            return context->throwError(QScriptContext::TypeError, QObject::tr("approx(): third argument is not a number"));
        eps = context->argument(2).toNumber();
        if (eps < 0)
            return context->throwError(QScriptContext::TypeError, QObject::tr("approx(): third argument must be greater than zero"));
    }
    if (!context->argument(0).isNumber())
        return context->throwError(QScriptContext::TypeError, QObject::tr("approx(): first argument is not a number"));
    if (!context->argument(1).isNumber())
        return context->throwError(QScriptContext::TypeError, QObject::tr("approx(): second argument is not a number"));
    double a = context->argument(0).toNumber();
    double b = context->argument(1).toNumber();

    return fabs(a - b) < eps;
}

QScriptValue createPoint2D(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() != 2)
        return context->throwError(QObject::tr("Point2D() takes exactly two arguments: Point2D(x, y)"));
    if (!context->argument(0).isNumber())
        return context->throwError(QScriptContext::TypeError, QObject::tr("Point2D(): first argument is not a number"));
    if (!context->argument(1).isNumber())
        return context->throwError(QScriptContext::TypeError, QObject::tr("Point2D(): second argument is not a number"));
    double x = context->argument(0).toNumber();
    double y = context->argument(1).toNumber();
    return engine->newQObject(new QPoint2D(x, y), QScriptEngine::ScriptOwnership);
}

QScriptValue toScriptValuePoint3D(QScriptEngine *engine, const msh::Point3D &point)
{
    QScriptValue obj = engine->newObject();
    obj.setProperty("x", point.x());
    obj.setProperty("y", point.y());
    obj.setProperty("z", point.z());
    return obj;
}

void fromScriptValuePoint3D(const QScriptValue &value, msh::Point3D &point)
{
    point.set(value.property("x").toNumber(), value.property("y").toNumber(), value.property("z").toNumber());
}

QScriptValue createPoint3D(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() != 3)
        return context->throwError(QObject::tr("Point3D() takes exactly three arguments: Point3D(x, y, z)"));
    if (!context->argument(0).isNumber())
        return context->throwError(QScriptContext::TypeError, QObject::tr("Point3D(): first argument is not a number"));
    if (!context->argument(1).isNumber())
        return context->throwError(QScriptContext::TypeError, QObject::tr("Point3D(): second argument is not a number"));
    if (!context->argument(2).isNumber())
        return context->throwError(QScriptContext::TypeError, QObject::tr("Point3D(): third argument is not a number"));
    double x = context->argument(0).toNumber();
    double y = context->argument(1).toNumber();
    double z = context->argument(2).toNumber();
    msh::Point3D point(x, y, z);
    return engine->toScriptValue(point);
}


QScriptValue printStd(QScriptContext *context, QScriptEngine *engine)
{
    for (int i = 0; i < context->argumentCount(); ++i) {
        if (i > 0)
            std::cout << ' ';
        std::cout << context->argument(i).toString().toAscii().data() ;
    }

    std::cout << std::endl;

    return engine->undefinedValue();
}

QScriptValue sum(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argument(0).isNumber())
        std::cout << "Number" << std::endl;
    if (context->argument(0).isObject())
    {
        QPoint2D *point = qscriptvalue_cast<QPoint2D*>(context->argument(0));
        if ( point != NULL )
            std::cout << "Point2D" << std::endl;
    }
    if (context->argument(0).isQMetaObject())
        std::cout << "QMetaObject" << std::endl;
    if (context->argument(0).isQObject())
        std::cout << "QObject" << std::endl;
}
