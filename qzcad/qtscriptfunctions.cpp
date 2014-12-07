#include <math.h>
#include <QObject>
#include <iostream>
#include "qtscriptfunctions.h"

QScriptValue approx(QScriptContext *ctx, QScriptEngine *eng)
{
    Q_UNUSED(eng);
    double eps = 1.0E-6;
    if (ctx->argumentCount() != 2 && ctx->argumentCount() != 3)
    {
        return ctx->throwError(QObject::tr("approx() takes two or three arguments: approx(x, y) or approx(x, y, eps)"));
    }
    if (ctx->argumentCount() == 3)
    {
        if (!ctx->argument(2).isNumber())
            return ctx->throwError(QScriptContext::TypeError, QObject::tr("approx(): third argument is not a number"));
        eps = ctx->argument(2).toNumber();
        if (eps < 0)
            return ctx->throwError(QScriptContext::TypeError, QObject::tr("approx(): third argument must be greater than zero"));
    }
    if (!ctx->argument(0).isNumber())
        return ctx->throwError(QScriptContext::TypeError, QObject::tr("approx(): first argument is not a number"));
    if (!ctx->argument(1).isNumber())
        return ctx->throwError(QScriptContext::TypeError, QObject::tr("approx(): second argument is not a number"));
    double a = ctx->argument(0).toNumber();
    double b = ctx->argument(1).toNumber();

    return fabs(a - b) < eps;
}

QScriptValue toScriptValuePoint2D(QScriptEngine *eng, const msh::Point2D &point)
{
    QScriptValue obj = eng->newObject();
    obj.setProperty("x", point.x());
    obj.setProperty("y", point.y());
    return obj;
}

void fromScriptValuePoint2D(const QScriptValue &value, msh::Point2D &point)
{
    point.set(value.property("x").toNumber(), value.property("y").toNumber());
}

QScriptValue createPoint2D(QScriptContext *ctx, QScriptEngine *eng)
{
    if (ctx->argumentCount() != 2)
        return ctx->throwError(QObject::tr("Point2D() takes exactly two arguments: Point2D(x, y)"));
    if (!ctx->argument(0).isNumber())
        return ctx->throwError(QScriptContext::TypeError, QObject::tr("Point2D(): first argument is not a number"));
    if (!ctx->argument(1).isNumber())
        return ctx->throwError(QScriptContext::TypeError, QObject::tr("Point2D(): second argument is not a number"));
    double x = ctx->argument(0).toNumber();
    double y = ctx->argument(1).toNumber();
    msh::Point2D point(x, y);
    return eng->toScriptValue(point);
}

QScriptValue toScriptValuePoint3D(QScriptEngine *eng, const msh::Point3D &point)
{
    QScriptValue obj = eng->newObject();
    obj.setProperty("x", point.x());
    obj.setProperty("y", point.y());
    obj.setProperty("z", point.z());
    return obj;
}

void fromScriptValuePoint3D(const QScriptValue &value, msh::Point3D &point)
{
    point.set(value.property("x").toNumber(), value.property("y").toNumber(), value.property("z").toNumber());
}

QScriptValue createPoint3D(QScriptContext *ctx, QScriptEngine *eng)
{
    if (ctx->argumentCount() != 3)
        return ctx->throwError(QObject::tr("Point3D() takes exactly three arguments: Point3D(x, y, z)"));
    if (!ctx->argument(0).isNumber())
        return ctx->throwError(QScriptContext::TypeError, QObject::tr("Point3D(): first argument is not a number"));
    if (!ctx->argument(1).isNumber())
        return ctx->throwError(QScriptContext::TypeError, QObject::tr("Point3D(): second argument is not a number"));
    if (!ctx->argument(2).isNumber())
        return ctx->throwError(QScriptContext::TypeError, QObject::tr("Point3D(): third argument is not a number"));
    double x = ctx->argument(0).toNumber();
    double y = ctx->argument(1).toNumber();
    double z = ctx->argument(2).toNumber();
    msh::Point3D point(x, y, z);
    return eng->toScriptValue(point);
}


QScriptValue printStd(QScriptContext *ctx, QScriptEngine *eng)
{
    for (int i = 0; i < ctx->argumentCount(); ++i) {
        if (i > 0)
            std::cout << ' ';
        std::cout << ctx->argument(i).toString().toAscii().data() ;
    }

    std::cout << std::endl;

    return eng->undefinedValue();
}
