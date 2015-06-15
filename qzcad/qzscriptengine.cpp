#include <iostream>
#include <math.h>
#include <list>
#include "qpoint2d.h"
#include "qpoint3d.h"
#include "qquadrilateralmesh2d.h"
#include "qtrianglemesh2d.h"

#include "qzscriptengine.h"

//using namespace msh;

double QZScriptEngine::epsilon_ = 1.0E-6;
Mesh *QZScriptEngine::mesh_ = NULL;

QZScriptEngine::QZScriptEngine(QObject *parent) :
    QScriptEngine(parent)
{
    mesh_ = NULL;
    // Операция приблизительного равентсва
    QScriptValue qsApprox = newFunction(approx);
    globalObject().setProperty("approx", qsApprox);
    // Переопределение печати
    QScriptValue qsPrint = newFunction(printStd);
    globalObject().setProperty("print", qsPrint);
    // Функция суммирования
    QScriptValue qsSum = newFunction(sum);
    globalObject().setProperty("sum", qsSum);
    // Функция конъюнкции
    QScriptValue qsCon = newFunction(con);
    globalObject().setProperty("con", qsCon);
    // Функция дизъюнкции
    QScriptValue qsDis = newFunction(dis);
    globalObject().setProperty("dis", qsDis);
    // Функция разности
    QScriptValue qsDiff = newFunction(diff);
    globalObject().setProperty("diff", qsDiff);
    // Точка на плоскости
    QScriptValue qsCreatePoint2D = newFunction(createPoint2D);
    globalObject().setProperty("Point2D", qsCreatePoint2D);
    // Точка в пространстве
    QScriptValue qsCreatePoint3D = newFunction(createPoint3D);
    globalObject().setProperty("Point3D", qsCreatePoint3D);
    // Двумерная сетка четырехугольников
    QScriptValue qsCreateQuadrilateralMesh2D = newFunction(createQuadrilateralMesh2D);
    globalObject().setProperty("Quads2D", qsCreateQuadrilateralMesh2D);
    // Двумерная сетка треугольников
    QScriptValue qsCreateTriangleMesh2D = newFunction(createTriangleMesh2D);
    globalObject().setProperty("Triangles2D", qsCreateTriangleMesh2D);
    // setMesh
    QScriptValue qsSetMesh = newFunction(setMesh);
    globalObject().setProperty("setMesh", qsSetMesh);
    // About
    QScriptValue qsAbout = newFunction(about);
    globalObject().setProperty("About", qsAbout);
}

double QZScriptEngine::epsilon() const
{
    return epsilon_;
}

void QZScriptEngine::setEpsilon(double epsilon)
{
    epsilon_ = epsilon;
}

Mesh *QZScriptEngine::mesh()
{
    return mesh_;
}

unsigned long QZScriptEngine::getNodeValuesSize() const
{
    return nodeValues_.size();
}

NamedFloatingVector &QZScriptEngine::getNodeValues(const unsigned long &i)
{
    return nodeValues_[i];
}

unsigned long QZScriptEngine::getElementValuesSize() const
{
    return elementValues_.size();
}

NamedFloatingVector &QZScriptEngine::getElementValues(const unsigned long &i)
{
    return elementValues_[i];
}

QScriptValue QZScriptEngine::about(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(context);
    QString message = tr("QZ Script Engine, version 1.0.0. Copyright 2014 Sergey Choporov. All rights reserved.");
    std::cout << message.toUtf8().data() << std::endl;
    return engine->undefinedValue();
}

QScriptValue QZScriptEngine::approx(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);
    double eps = epsilon_;
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

QScriptValue QZScriptEngine::createPoint2D(QScriptContext *context, QScriptEngine *engine)
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

QScriptValue QZScriptEngine::createPoint3D(QScriptContext *context, QScriptEngine *engine)
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
    return engine->newQObject(new QPoint3D(x, y, z), QScriptEngine::ScriptOwnership);
}

QScriptValue QZScriptEngine::createQuadrilateralMesh2D(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() == 4)
    {
        QString typeError = QObject::tr("Quads2D(count: Integer, center: Point2D, radius: Floating, part: {1, 2, 4}): argument type error (%1).");
        if (!context->argument(0).isNumber())
            return context->throwError(typeError.arg("count"));
        if (!context->argument(1).isQObject() || qscriptvalue_cast<QPoint2D *>(context->argument(1)) == NULL)
            return context->throwError(typeError.arg("center"));
        if (!context->argument(2).isNumber())
            return context->throwError(typeError.arg("radius"));
        if (!context->argument(3).isNumber())
            return context->throwError(typeError.arg("part"));
        UInteger count = context->argument(0).toUInt32();
        QPoint2D *center = qscriptvalue_cast<QPoint2D *>(context->argument(1));
        double radius = context->argument(2).toNumber();
        int part = context->argument(3).toNumber();
        return engine->newQObject(new QQuadrilateralMesh2D(count, msh::Point2D(center->x(), center->y()), radius, part), QScriptEngine::ScriptOwnership);
    }
    if (context->argumentCount() == 5)
    {
        QString typeError = QObject::tr("Quads2D(xCount: Integer, yCount: Integer, origin: Point2D, width: Floating, height: Floating): argument type error (%1).");
        if (!context->argument(0).isNumber())
            return context->throwError(typeError.arg("xCount"));
        if (!context->argument(1).isNumber())
            return context->throwError(typeError.arg("yCount"));
        if (!context->argument(2).isQObject() || qscriptvalue_cast<QPoint2D *>(context->argument(2)) == NULL)
            return context->throwError(typeError.arg("origin"));
        if (!context->argument(3).isNumber())
            return context->throwError(typeError.arg("width"));
        if (!context->argument(4).isNumber())
            return context->throwError(typeError.arg("height"));
        UInteger xCount = context->argument(0).toUInt32();
        UInteger yCount = context->argument(1).toUInt32();
        QPoint2D *origin = qscriptvalue_cast<QPoint2D *>(context->argument(2));
        double width = context->argument(3).toNumber();
        double height = context->argument(4).toNumber();
        return engine->newQObject(new QQuadrilateralMesh2D(xCount, yCount, origin->x(), origin->y(), width, height), QScriptEngine::ScriptOwnership);
    }
    return context->throwError(QObject::tr("Quads2D(): arguments count error."));
}

QScriptValue QZScriptEngine::createTriangleMesh2D(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() == 5)
    {
        QString typeError = QObject::tr("Triangles2D(xCount: Integer, yCount: Integer, origin: Point2D, width: Floating, height: Floating): argument type error (%1).");
        if (!context->argument(0).isNumber())
            return context->throwError(typeError.arg("xCount"));
        if (!context->argument(1).isNumber())
            return context->throwError(typeError.arg("yCount"));
        if (!context->argument(2).isQObject() || qscriptvalue_cast<QPoint2D *>(context->argument(2)) == NULL)
            return context->throwError(typeError.arg("origin"));
        if (!context->argument(3).isNumber())
            return context->throwError(typeError.arg("width"));
        if (!context->argument(4).isNumber())
            return context->throwError(typeError.arg("height"));
        UInteger xCount = context->argument(0).toUInt32();
        UInteger yCount = context->argument(1).toUInt32();
        QPoint2D *origin = qscriptvalue_cast<QPoint2D *>(context->argument(2));
        double width = context->argument(3).toNumber();
        double height = context->argument(4).toNumber();
        return engine->newQObject(new QTriangleMesh2D(xCount, yCount, origin->x(), origin->y(), width, height), QScriptEngine::ScriptOwnership);
    }
    else if (context->argumentCount() == 6 || context->argumentCount() == 7)
    {
        QString typeError = QObject::tr("Triangles2D(xCount: Integer, yCount: Integer, origin: Point2D, width: Floating, height: Floating, function: Function[, points: Array]): argument type error (%1).");
        std::list<msh::Point2D> pointList;
        if (!context->argument(0).isNumber())
            return context->throwError(typeError.arg("xCount"));
        if (!context->argument(1).isNumber())
            return context->throwError(typeError.arg("yCount"));
        if (!context->argument(2).isQObject() || qscriptvalue_cast<QPoint2D *>(context->argument(2)) == NULL)
            return context->throwError(typeError.arg("origin"));
        if (!context->argument(3).isNumber())
            return context->throwError(typeError.arg("width"));
        if (!context->argument(4).isNumber())
            return context->throwError(typeError.arg("height"));
        QScriptValue function = context->argument(5);
        if (!function.isFunction())
            return context->throwError(typeError.arg("function"));

        UInteger xCount = context->argument(0).toUInt32();
        UInteger yCount = context->argument(1).toUInt32();
        QPoint2D *origin = qscriptvalue_cast<QPoint2D *>(context->argument(2));
        double width = context->argument(3).toNumber();
        double height = context->argument(4).toNumber();
        // функция для вычисления квадрата числа (C++0x)
        auto func = [&](double x, double y)
        {
            QScriptValueList args;
            args << x << y;
            return function.call(QScriptValue(), args).toNumber();
        };
        if (context->argumentCount() == 7)
        {
            if (!context->argument(6).isArray())
                return context->throwError(typeError.arg("points"));
            QScriptValue array = context->argument(6);
            for (int i = 0; i < array.property("length").toInteger(); i++)
            {
                QPoint2D *point = qscriptvalue_cast<QPoint2D *>(array.property(i));
                pointList.push_back(Point2D(point->x(), point->y()));
            }
        }
        return engine->newQObject(new QTriangleMesh2D(xCount, yCount, origin->x(), origin->y(), width, height, func, pointList), QScriptEngine::ScriptOwnership);
    }
    return context->throwError(QObject::tr("Triangles2D(xCount: Integer, yCount: Integer, origin: Point2D, width: Floating, height: Floating): arguments count error."));
}

QScriptValue QZScriptEngine::printStd(QScriptContext *context, QScriptEngine *engine)
{
    for (register int i = 0; i < context->argumentCount(); ++i)
    {
        if (i > 0)
            std::cout << ' ';
        std::cout << context->argument(i).toString().toStdString() ;
    }

    std::cout << std::endl;

    return engine->undefinedValue();
}

QScriptValue QZScriptEngine::setMesh(QScriptContext *context, QScriptEngine *engine)
{
    QString typeError = QObject::tr("setMesh(mesh): argument must be a mesh object");
    if (mesh_ != NULL)
    {
        delete mesh_;
        mesh_ = NULL;
    }
    if (context->argumentCount() != 1)
        return context->throwError(QObject::tr("setMesh() takes exactly one argument: setMesh(mesh)"));
    if (!context->argument(0).isQObject())
        return context->throwError(typeError);
    if (qscriptvalue_cast<QQuadrilateralMesh2D *>(context->argument(0)) != NULL)
    {
        mesh_ = new QuadrilateralMesh2D(qscriptvalue_cast<QQuadrilateralMesh2D *>(context->argument(0)));
    }
    else if (qscriptvalue_cast<QTriangleMesh2D *>(context->argument(0)) != NULL)
    {
        mesh_ = new TriangleMesh2D(qscriptvalue_cast<QTriangleMesh2D *>(context->argument(0)));
    }
    else
    {
        return context->throwError(typeError);
    }
    return engine->undefinedValue();
}

QScriptValue QZScriptEngine::sum(QScriptContext *context, QScriptEngine *engine)
{
    QString typeError = QObject::tr("sum(a, b, c, ...): all arguments must be a same type: %1");
    if (context->argument(0).isNumber())
    {
        for (int i = 1; i < context->argumentCount(); i++)
            if (!context->argument(i).isNumber())
                return context->throwError(typeError.arg("number"));
        double sum = context->argument(0).toNumber();
        for (int i = 1; i < context->argumentCount(); i++)
            sum += context->argument(i).toNumber();
        return sum;
    }
    if (context->argument(0).isQObject())
    {
        if ( qscriptvalue_cast<QPoint2D*>(context->argument(0)) != NULL )
        {
            for (int i = 1; i < context->argumentCount(); i++)
                if (qscriptvalue_cast<QPoint2D*>(context->argument(i)) == NULL)
                    return context->throwError(typeError.arg("Point2D"));
            QPoint2D *p0 = qscriptvalue_cast<QPoint2D*>(context->argument(0));
            QPoint2D *sum = new QPoint2D(*p0);
            for (int i = 1; i < context->argumentCount(); i++)
            {
                QPoint2D *pi = qscriptvalue_cast<QPoint2D*>(context->argument(i));
                sum->setX(sum->x() + pi->x());
                sum->setY(sum->y() + pi->y());
            }
            return engine->newQObject(sum, QScriptEngine::ScriptOwnership);
        }
        if ( qscriptvalue_cast<QPoint3D*>(context->argument(0)) != NULL )
        {
            for (int i = 1; i < context->argumentCount(); i++)
                if (qscriptvalue_cast<QPoint3D*>(context->argument(i)) == NULL)
                    return context->throwError(typeError.arg("Point3D"));
            QPoint3D *p0 = qscriptvalue_cast<QPoint3D*>(context->argument(0));
            QPoint3D *sum = new QPoint3D(*p0);
            for (int i = 1; i < context->argumentCount(); i++)
            {
                QPoint3D *pi = qscriptvalue_cast<QPoint3D*>(context->argument(i));
                sum->setX(sum->x() + pi->x());
                sum->setY(sum->y() + pi->y());
                sum->setZ(sum->z() + pi->z());
            }
            return engine->newQObject(sum, QScriptEngine::ScriptOwnership);
        }
    }
    return engine->undefinedValue();
}

QScriptValue QZScriptEngine::con(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);
    QString typeError = QObject::tr("con(a, b, c, ...): all arguments must have type Number. Argument # %1.");
    for (int i = 0; i < context->argumentCount(); i++)
        if (!context->argument(i).isNumber())
            return context->throwError(typeError.arg(i + 1));
    double con = context->argument(0).toNumber();
    for (int i = 1; i < context->argumentCount(); i++)
    {
        double x = context->argument(i).toNumber();
        con += (x - sqrt(con*con + x*x));
    }
    return con;
}

QScriptValue QZScriptEngine::dis(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);
    QString typeError = QObject::tr("dis(a, b, c, ...): all arguments must have type Number. Argument # %1.");
    for (int i = 0; i < context->argumentCount(); i++)
        if (!context->argument(i).isNumber())
            return context->throwError(typeError.arg(i + 1));
    double dis = context->argument(0).toNumber();
    for (int i = 1; i < context->argumentCount(); i++)
    {
        double x = context->argument(i).toNumber();
        dis += (x + sqrt(dis*dis + x*x));
    }
    return dis;
}

QScriptValue QZScriptEngine::diff(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);
    QString typeError = QObject::tr("diff(a, b, c, ...): all arguments must have type Number. Argument # %1.");
    for (int i = 0; i < context->argumentCount(); i++)
        if (!context->argument(i).isNumber())
            return context->throwError(typeError.arg(i + 1));
    double diff = context->argument(0).toNumber();
    for (int i = 1; i < context->argumentCount(); i++)
    {
        double x = -context->argument(i).toNumber();
        diff += (x - sqrt(diff*diff + x*x));
    }
    return diff;
}


