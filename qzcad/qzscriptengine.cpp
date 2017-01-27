#include <iostream>
#include <math.h>
#include <list>
#include "qpoint2d.h"
#include "qpoint3d.h"
#include "qquadrilateralmesh2d.h"
#include "qtrianglemesh2d.h"
#include "qquadrilateralmesh3d.h"
#include "qsegmentmesh2d.h"
#include "qtrianglemesh3d.h"

#include "rfunctions.h"

#include "qfemcondition.h"

#include "planestressstrain.h"
#include "mindlinplatebending.h"
#include "mindlinplatelaminated.h"
#include "mindlinshellbending.h"
#include "mindlinshelllaminated.h"

#include "qzscriptengine.h"

//using namespace msh;

double QZScriptEngine::epsilon_ = 1.0E-6;
Mesh *QZScriptEngine::mesh_ = NULL;
Fem *QZScriptEngine::fem_ = NULL;

QZScriptEngine::QZScriptEngine(QObject *parent) :
    QScriptEngine(parent)
{
    mesh_ = NULL;
    fem_ = NULL;
    // Операция приблизительного равентсва
    globalObject().setProperty("approx", newFunction(approx));
    // Переопределение печати
    globalObject().setProperty("print", newFunction(printStd));
    // Функция суммирования
    globalObject().setProperty("sum", newFunction(sum));
    // Функция конъюнкции
    globalObject().setProperty("con", newFunction(con));
    // Функция дизъюнкции
    globalObject().setProperty("dis", newFunction(dis));
    // Функция разности
    globalObject().setProperty("diff", newFunction(diff));
    // circle
    globalObject().setProperty("circle", newFunction(circle));
    // ellipse
    globalObject().setProperty("ellipse", newFunction(ellipse));
    // band
    globalObject().setProperty("band", newFunction(band));
    // line
    globalObject().setProperty("line", newFunction(line));
    // rectangle
    globalObject().setProperty("rectangle", newFunction(rectangle));
    // convex
    globalObject().setProperty("convex", newFunction(convex));
    // regular
    globalObject().setProperty("regular", newFunction(regular));
    // ellipsoid
    globalObject().setProperty("ellipsoid", newFunction(ellipsoid));
    // sphere
    globalObject().setProperty("sphere", newFunction(sphere));
    // plane
    globalObject().setProperty("plane", newFunction(plane));
    // cuboid
    globalObject().setProperty("cuboid", newFunction(cuboid));

    // Точка (плоскость или пространство)
    globalObject().setProperty("Point", newFunction(createPoint));

    // Двумерная сетка четырехугольников
    globalObject().setProperty("Quads2D", newFunction(createQuadrilateralMesh2D));
    // Двумерная сетка отрезков
    globalObject().setProperty("Segments2D", newFunction(createSegmentMesh2D));
    // Двумерная картина линий уровня
    globalObject().setProperty("ContourGraph", newFunction(createContourGraph));
    // Двумерная сетка треугольников
    globalObject().setProperty("Triangles2D", newFunction(createTriangleMesh2D));
    // Триангуляция Делоне
    globalObject().setProperty("delaunay", newFunction(createDelaunay));
    // Триангуляция Делоне с использованием сглаживания методом Рапперта
    globalObject().setProperty("ruppert", newFunction(createRuppert));
    // Поверхностная сетка четырехугольников
    globalObject().setProperty("ShellQuads", newFunction(createQuadrilateralMesh3D));
    // Поверхностная сетка четырехугольников в цилиндрических координатах
    globalObject().setProperty("CylinderQuads", newFunction(createCylinderQuads));
    // Поверхностная сетка четырехугольников в конических координатах
    globalObject().setProperty("ConeQuads", newFunction(createConeQuads));
    // Поверхностная сетка треугольников
    globalObject().setProperty("ShellTriangles", newFunction(createTriangleMesh3D));
    // Поверхностная сетка треугольников
    globalObject().setProperty("CylinderTriangles", newFunction(createCylinderTriangles));
    // Поверхностная сетка треугольников в конических координатах
    globalObject().setProperty("ConeTriangles", newFunction(createConeTriangles));
    // Поверхностная сетка треугольников в параметрических координатах
    globalObject().setProperty("ParametricTriangles", newFunction(createParametricTriangles));
    // Марширующие кубики
    globalObject().setProperty("MarchingCubes", newFunction(createMarchingCubes));

    // setMesh
    globalObject().setProperty("setMesh", newFunction(setMesh));
    // currentMesh
    globalObject().setProperty("currentMesh", newFunction(currentMesh));

    // About
    globalObject().setProperty("About", newFunction(about));

    /// FEM section
    globalObject().setProperty("PlaneStress", newFunction(planeStress));
    globalObject().setProperty("PlaneStrain", newFunction(planeStrain));

    globalObject().setProperty("BoundaryCondition", newFunction(createBoundaryCondition));
    globalObject().setProperty("NodalForce", newFunction(createNodalForce));
    QScriptValue qsSurfaceForce = newFunction(createSurfaceForce);
    globalObject().setProperty("SurfaceForce", qsSurfaceForce);
    globalObject().setProperty("EdgeForce", qsSurfaceForce); // создание синонима
    QScriptValue qsVolumeForce = newFunction(createVolumeForce);
    globalObject().setProperty("VolumeForce", qsVolumeForce);
    globalObject().setProperty("PlateDistributedForce", qsVolumeForce); // создание синонима
    globalObject().setProperty("ShellDistributedForce", qsVolumeForce); // создание синонима
    // направления действия нагрузки
    globalObject().setProperty("ALL", FemCondition::ALL);
    globalObject().setProperty("FIRST", FemCondition::FIRST);
    globalObject().setProperty("SECOND", FemCondition::SECOND);
    globalObject().setProperty("THIRD", FemCondition::THIRD);
    globalObject().setProperty("FOURTH", FemCondition::FOURTH);
    globalObject().setProperty("FIFTH", FemCondition::FIFTH);
    globalObject().setProperty("SIXTH", FemCondition::SIXTH);
    globalObject().setProperty("SEVENTH", FemCondition::SEVENTH);
    globalObject().setProperty("EIGHTH", FemCondition::EIGHTH);
    globalObject().setProperty("NINETH", FemCondition::NINETH);
    globalObject().setProperty("TENTH", FemCondition::TENTH);

    globalObject().setProperty("MindlinPlate", newFunction(mindlinPlate));

    globalObject().setProperty("MindlinShell", newFunction(mindlinShell));

    globalObject().setProperty("values", newFunction(reportValues));
}

double QZScriptEngine::epsilon()
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

void QZScriptEngine::setMesh(Mesh *mesh)
{
    if (mesh_ != NULL)
        delete mesh_;

    mesh_ = mesh;
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

QScriptValue QZScriptEngine::createPoint(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() == 2)
    {
        if (!context->argument(0).isNumber())
            return context->throwError(QScriptContext::TypeError, QObject::tr("Point(): first argument is not a number"));
        if (!context->argument(1).isNumber())
            return context->throwError(QScriptContext::TypeError, QObject::tr("Point(): second argument is not a number"));
        double x = context->argument(0).toNumber();
        double y = context->argument(1).toNumber();
        return engine->newQObject(new QPoint2D(x, y), QScriptEngine::ScriptOwnership);
    }
    else if (context->argumentCount() == 3)
    {
        if (!context->argument(0).isNumber())
            return context->throwError(QScriptContext::TypeError, QObject::tr("Point(): first argument is not a number"));
        if (!context->argument(1).isNumber())
            return context->throwError(QScriptContext::TypeError, QObject::tr("Point(): second argument is not a number"));
        if (!context->argument(2).isNumber())
            return context->throwError(QScriptContext::TypeError, QObject::tr("Point(): third argument is not a number"));
        double x = context->argument(0).toNumber();
        double y = context->argument(1).toNumber();
        double z = context->argument(2).toNumber();
        return engine->newQObject(new QPoint3D(x, y, z), QScriptEngine::ScriptOwnership);
    }
    return context->throwError(QObject::tr("Point() takes two or three arguments: Point(x, y) or Point(x, y, z)"));
}

QScriptValue QZScriptEngine::createQuadrilateralMesh2D(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() == 4)
    {
        QString typeError = QObject::tr("Quads2D(count: Integer, center: Point, radius: Floating, part: {1, 2, 4}): argument type error (%1).");
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

        QQuadrilateralMesh2D *qmo = new QQuadrilateralMesh2D();
        qmo->circleDomain(count, msh::Point2D(center->x(), center->y()), radius, part);

        return engine->newQObject(qmo, QScriptEngine::ScriptOwnership);
    }
    if (context->argumentCount() == 5)
    {
        QString typeError = QObject::tr("Quads2D(xCount: Integer, yCount: Integer, origin: Point, width: Floating, height: Floating): argument type error (%1).");
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

        QQuadrilateralMesh2D *qmo = new QQuadrilateralMesh2D();
        qmo->rectangleDomain(xCount, yCount, origin->x(), origin->y(), width, height);

        return engine->newQObject(qmo, QScriptEngine::ScriptOwnership);
    }

    if (context->argumentCount() == 6 || context->argumentCount() == 7)
    {
        QString typeError = QObject::tr("Quads2D(xCount: Integer, yCount: Integer, origin: Point, width: Floating, height: Floating, function: Function[, points: Array]): argument type error (%1).");
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

        QQuadrilateralMesh2D *qmo = new QQuadrilateralMesh2D();
        qmo->functionalDomain(xCount, yCount, origin->x(), origin->y(), width, height, func, pointList);

        return engine->newQObject(qmo, QScriptEngine::ScriptOwnership);
    }
    return context->throwError(QObject::tr("Quads2D(): arguments count error."));
}

QScriptValue QZScriptEngine::createSegmentMesh2D(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() == 4)
    {
        QString typeError = QObject::tr("Segments2D(count: Integer, min: Floating, max: Floating, function: Function): argument type error (%1).");
        if (!context->argument(0).isNumber())
            return context->throwError(typeError.arg("count"));
        if (!context->argument(1).isNumber())
            return context->throwError(typeError.arg("min"));
        if (!context->argument(2).isNumber())
            return context->throwError(typeError.arg("max"));
        QScriptValue function = context->argument(3);
        if (!function.isFunction())
            return context->throwError(typeError.arg("function"));
        auto domain = [&](double t)
        {
            QScriptValueList args;
            args << t;
            QPoint2D *p = qscriptvalue_cast<QPoint2D *>(function.call(QScriptValue(), args));
            return Point2D(*p);
        };
        UInteger count = context->argument(0).toUInt32();
        double tmin = context->argument(1).toNumber();
        double tmax = context->argument(2).toNumber();
        QSegmentMesh2D *smo = new QSegmentMesh2D();
        smo->parametricDomain(count, tmin, tmax, domain);
        return engine->newQObject(smo, QScriptEngine::ScriptOwnership);
    }
    else if (context->argumentCount() == 6 || (context->argumentCount() >= 7 && context->argumentCount() <= 9))
    {
        QString typeError = QObject::tr("Segments2D(xCount: Integer, yCount: Integer, origin: Point, width: Floating, height: Floating, function: Function[, points: Array]): argument type error (%1).");
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

        if (context->argumentCount() == 7 && context->argument(6).isArray())
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
        else if (context->argumentCount() >= 7)
        {
            // случай контакта
            double delta = -1.0; // парметр окрестности сгущения: по умолчани не сгущаем
            QString typeError = QObject::tr("Segments2D(xCount: Integer, yCount: Integer, origin: Point, width: Floating, height: Floating, func_a: Function, func_b: Function[, points: Array, delta: Floating]): argument type error (%1).");
            if (!context->argument(6).isFunction())
                return context->throwError(typeError.arg("func_b"));
            QScriptValue function_b = context->argument(6);
            auto func_b = [&](double x, double y)
            {
                QScriptValueList args;
                args << x << y;
                return function_b.call(QScriptValue(), args).toNumber();
            };
            if (context->argumentCount() >= 8)
            {
                if (!context->argument(7).isArray())
                    return context->throwError(typeError.arg("points"));
                QScriptValue array = context->argument(7);
                for (int i = 0; i < array.property("length").toInteger(); i++)
                {
                    QPoint2D *point = qscriptvalue_cast<QPoint2D *>(array.property(i));
                    pointList.push_back(Point2D(point->x(), point->y()));
                }
            }
            if (context->argumentCount() == 9)
            {
                if (!context->argument(8).isNumber())
                    return context->throwError(typeError.arg("delta"));
                delta = context->argument(8).toNumber();
                std::cout << delta << std::endl;
            }
            QSegmentMesh2D *smo = new QSegmentMesh2D();
            smo->functionalDomain(xCount, yCount, origin->x(), origin->y(), width, height, func, func_b, pointList, delta);

            return engine->newQObject(smo, QScriptEngine::ScriptOwnership);
        }

        QSegmentMesh2D *smo = new QSegmentMesh2D();
        smo->functionalDomain(xCount, yCount, origin->x(), origin->y(), width, height, func, pointList);

        return engine->newQObject(smo, QScriptEngine::ScriptOwnership);
    }
    return context->throwError(QObject::tr("Segments2D(xCount: Integer, yCount: Integer, origin: Point, width: Floating, height: Floating): arguments count error."));
}

QScriptValue QZScriptEngine::createContourGraph(QScriptContext *context, QScriptEngine *engine)
{
    if (7 <= context->argumentCount() && context->argumentCount() <= 8)
        {
            QString typeError = QObject::tr("ContourGraph(xCount: Integer, yCount: Integer, origin: Point, width: Floating, height: Floating, function: Function, contours: Integer[, points: Array]): argument type error (%1).");
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

            if (!context->argument(6).isNumber())
                return context->throwError(typeError.arg("contours"));
            int contours = context->argument(6).toInt32();

            if (context->argumentCount() == 8)
            {
                if (!context->argument(7).isArray())
                    return context->throwError(typeError.arg("points"));
                QScriptValue array = context->argument(7);
                for (int i = 0; i < array.property("length").toInteger(); i++)
                {
                    QPoint2D *point = qscriptvalue_cast<QPoint2D *>(array.property(i));
                    pointList.push_back(Point2D(point->x(), point->y()));
                }
            }

            QSegmentMesh2D *smo = new QSegmentMesh2D();
            smo->contourGraph(xCount, yCount, origin->x(), origin->y(), width, height, func, pointList, contours, true);

            return engine->newQObject(smo, QScriptEngine::ScriptOwnership);
        }
    return context->throwError(QObject::tr("ContourGraph(xCount: Integer, yCount: Integer, origin: Point, width: Floating, height: Floating, function: Function, contours: Integer[, points: Array]): arguments count error."));
}

QScriptValue QZScriptEngine::createTriangleMesh2D(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() == 1)
    {
        QString typeError = QObject::tr("Triangles2D(mesh: Mesh2D): argument type error.");
        if (!context->argument(0).isQObject())
            return context->throwError(typeError);
        if (qscriptvalue_cast<QSegmentMesh2D *>(context->argument(0)) != NULL)
        {
            SegmentMesh2D sm(qscriptvalue_cast<QSegmentMesh2D *>(context->argument(0)));

            QTriangleMesh2D *tmo = new QTriangleMesh2D();
            tmo->delaunay(&sm);

            return engine->newQObject(tmo, QScriptEngine::ScriptOwnership);
        }
    }
    else if (context->argumentCount() == 5)
    {
        QString typeError = QObject::tr("Triangles2D(xCount: Integer, yCount: Integer, origin: Point, width: Floating, height: Floating): argument type error (%1).");
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

        QTriangleMesh2D *tmo = new QTriangleMesh2D();
        tmo->rectangleDomain(xCount, yCount, origin->x(), origin->y(), width, height);

        return engine->newQObject(tmo, QScriptEngine::ScriptOwnership);
    }
    else if (context->argumentCount() == 6 || context->argumentCount() == 7)
    {
        QString typeError = QObject::tr("Triangles2D(xCount: Integer, yCount: Integer, origin: Point, width: Floating, height: Floating, function: Function[, points: Array]): argument type error (%1).");
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

        QTriangleMesh2D *tmo = new QTriangleMesh2D();
        tmo->functionalDomain(xCount, yCount, origin->x(), origin->y(), width, height, func, pointList);

        return engine->newQObject(tmo, QScriptEngine::ScriptOwnership);
    }
    return context->throwError(QObject::tr("Triangles2D(xCount: Integer, yCount: Integer, origin: Point, width: Floating, height: Floating): arguments count error."));
}

QScriptValue QZScriptEngine::createDelaunay(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() == 6 || context->argumentCount() == 7)
        {
            QString typeError = QObject::tr("delaunay(xCount: Integer, yCount: Integer, origin: Point, width: Floating, height: Floating, function: Function[, points: Array]): argument type error (%1).");
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
            QTriangleMesh2D  *qtm = new QTriangleMesh2D();
            qtm->delaunay(xCount, yCount, origin->x(), origin->y(), width, height, func, pointList);
            return engine->newQObject(qtm, QScriptEngine::ScriptOwnership);
        }
    return context->throwError(QObject::tr("delaunay(xCount: Integer, yCount: Integer, origin: Point, width: Floating, height: Floating): arguments count error."));
}

QScriptValue QZScriptEngine::createRuppert(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() == 6  || (context->argumentCount() >= 7 && context->argumentCount() <= 9))
        {
            QString typeError = QObject::tr("ruppert(xCount: Integer, yCount: Integer, origin: Point, width: Floating, height: Floating, function: Function[, points: Array]): argument type error (%1).");
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
            if (context->argumentCount() == 7 && context->argument(6).isArray())
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
            else if (context->argumentCount() >= 7)
            {
                // случай контакта
                double delta = -1.0; // парметр окрестности сгущения: по умолчани не сгущаем
                QString typeError = QObject::tr("Segments2D(xCount: Integer, yCount: Integer, origin: Point, width: Floating, height: Floating, func_a: Function, func_b: Function[, points: Array, delta: Floating]): argument type error (%1).");
                if (!context->argument(6).isFunction())
                    return context->throwError(typeError.arg("func_b"));
                QScriptValue function_b = context->argument(6);
                auto func_b = [&](double x, double y)
                {
                    QScriptValueList args;
                    args << x << y;
                    return function_b.call(QScriptValue(), args).toNumber();
                };
                if (context->argumentCount() >= 8)
                {
                    if (!context->argument(7).isArray())
                        return context->throwError(typeError.arg("points"));
                    QScriptValue array = context->argument(7);
                    for (int i = 0; i < array.property("length").toInteger(); i++)
                    {
                        QPoint2D *point = qscriptvalue_cast<QPoint2D *>(array.property(i));
                        pointList.push_back(Point2D(point->x(), point->y()));
                    }
                }
                if (context->argumentCount() == 9)
                {
                    if (!context->argument(8).isNumber())
                        return context->throwError(typeError.arg("delta"));
                    delta = context->argument(8).toNumber();
                    std::cout << delta << std::endl;
                }
                QTriangleMesh2D  *qtm = new QTriangleMesh2D();
                qtm->ruppert(xCount, yCount, origin->x(), origin->y(), width, height, func, func_b, pointList, delta);

                return engine->newQObject(qtm, QScriptEngine::ScriptOwnership);
            }
            QTriangleMesh2D  *qtm = new QTriangleMesh2D();
            qtm->ruppert(xCount, yCount, origin->x(), origin->y(), width, height, func, pointList);
            return engine->newQObject(qtm, QScriptEngine::ScriptOwnership);
        }
    return context->throwError(QObject::tr("ruppert(xCount: Integer, yCount: Integer, origin: Point, width: Floating, height: Floating): arguments count error."));
}

QScriptValue QZScriptEngine::createQuadrilateralMesh3D(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() == 4)
    {
        QString typeError = QObject::tr("ShellQuads(rCount: Integer, lCount: Integer, radius: Floating, length: Floating): argument type error (%1).");
        if (!context->argument(0).isNumber())
            return context->throwError(typeError.arg("rCount"));
        if (!context->argument(1).isNumber())
            return context->throwError(typeError.arg("lCount"));
        if (!(context->argument(2).isNumber() || context->argument(2).isFunction()))
            return context->throwError(typeError.arg("radius"));
        if (!context->argument(3).isNumber())
            return context->throwError(typeError.arg("length"));
        UInteger rCount = context->argument(0).toUInt32();
        UInteger lCount = context->argument(1).toUInt32();
        double length = context->argument(3).toNumber();

        QQuadrilateralMesh3D *qmo = new QQuadrilateralMesh3D();

        if (context->argument(2).isNumber())
        {
            double radius = context->argument(2).toNumber();
            qmo->cylinderDomain(rCount, lCount, radius, length);
        }
        else
        {
            QScriptValue radius = context->argument(2);
            auto radius_func = [&](double x)
            {
                QScriptValueList args;
                args << x;
                return radius.call(QScriptValue(), args).toNumber();
            };
            qmo->cylinderDomain(rCount, lCount, radius_func, length);
        }

        return engine->newQObject(qmo, QScriptEngine::ScriptOwnership);
    }
    return context->throwError(QObject::tr("ShellQuads(): arguments count error."));
}

QScriptValue QZScriptEngine::createCylinderQuads(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() == 4)
    {
        QString typeError = QObject::tr("CylinderQuads(rCount: Integer, lCount: Integer, radius: Floating, length: Floating): argument type error (%1).");
        if (!context->argument(0).isNumber())
            return context->throwError(typeError.arg("rCount"));
        if (!context->argument(1).isNumber())
            return context->throwError(typeError.arg("lCount"));
        if (!(context->argument(2).isNumber() || context->argument(2).isFunction()))
            return context->throwError(typeError.arg("radius"));
        if (!context->argument(3).isNumber())
            return context->throwError(typeError.arg("length"));
        UInteger rCount = context->argument(0).toUInt32();
        UInteger lCount = context->argument(1).toUInt32();
        double length = context->argument(3).toNumber();

        QQuadrilateralMesh3D *qmo = new QQuadrilateralMesh3D();

        if (context->argument(2).isNumber())
        {
            double radius = context->argument(2).toNumber();
            qmo->cylinderDomain(rCount, lCount, radius, length);
        }
        else
        {
            QScriptValue radius = context->argument(2);
            auto radius_func = [&](double x)
            {
                QScriptValueList args;
                args << x;
                return radius.call(QScriptValue(), args).toNumber();
            };
            qmo->cylinderDomain(rCount, lCount, radius_func, length);
        }

        return engine->newQObject(qmo, QScriptEngine::ScriptOwnership);
    }
    return context->throwError(QObject::tr("CylinderQuads(): arguments count error."));
}

QScriptValue QZScriptEngine::createConeQuads(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() == 5)
    {
        QString typeError = QObject::tr("ConeQuads(rCount: Integer, lCount: Integer, radiusBottom: Floating, radiusTop: Floating, length: Floating): argument type error (%1).");
        if (!context->argument(0).isNumber())
            return context->throwError(typeError.arg("rCount"));
        if (!context->argument(1).isNumber())
            return context->throwError(typeError.arg("lCount"));
        if (!context->argument(2).isNumber())
            return context->throwError(typeError.arg("radiusBottom"));
        if (!context->argument(3).isNumber())
            return context->throwError(typeError.arg("radiusTop"));
        if (!context->argument(4).isNumber())
            return context->throwError(typeError.arg("length"));
        UInteger rCount = context->argument(0).toUInt32();
        UInteger lCount = context->argument(1).toUInt32();
        double radiusBottom = context->argument(2).toNumber();
        double radiusTop = context->argument(3).toNumber();
        double length = context->argument(4).toNumber();

        QQuadrilateralMesh3D *qmo = new QQuadrilateralMesh3D();
        qmo->coneDomain(rCount, lCount, radiusBottom, radiusTop, length);

        return engine->newQObject(qmo, QScriptEngine::ScriptOwnership);
    }
    return context->throwError(QObject::tr("ConeQuads(): arguments count error."));
}

QScriptValue QZScriptEngine::createTriangleMesh3D(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() == 4)
    {
        QString typeError = QObject::tr("ShellTriangles(rCount: Integer, lCount: Integer, radius: Floating, length: Floating): argument type error (%1).");
        if (!context->argument(0).isNumber())
            return context->throwError(typeError.arg("rCount"));
        if (!context->argument(1).isNumber())
            return context->throwError(typeError.arg("lCount"));
        if (!context->argument(2).isNumber())
            return context->throwError(typeError.arg("radius"));
        if (!context->argument(3).isNumber())
            return context->throwError(typeError.arg("length"));
        UInteger rCount = context->argument(0).toUInt32();
        UInteger lCount = context->argument(1).toUInt32();
        double radius = context->argument(2).toNumber();
        double length = context->argument(3).toNumber();

        QTriangleMesh3D *tmo = new QTriangleMesh3D();
        tmo->cylinderDomain(rCount, lCount, radius, length);

        return engine->newQObject(tmo, QScriptEngine::ScriptOwnership);
    }
    return context->throwError(QObject::tr("ShellTriangles(): arguments count error."));
}

QScriptValue QZScriptEngine::createCylinderTriangles(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() == 4 || context->argumentCount() == 5)
    {
        QString typeError = QObject::tr("CylinderTriangles(rCount: Integer, lCount: Integer, radius: Floating, length: Floating[, func: Function]): argument type error (%1).");
        if (!context->argument(0).isNumber())
            return context->throwError(typeError.arg("rCount"));
        if (!context->argument(1).isNumber())
            return context->throwError(typeError.arg("lCount"));
        if (!context->argument(2).isNumber())
            return context->throwError(typeError.arg("radius"));
        if (!context->argument(3).isNumber())
            return context->throwError(typeError.arg("length"));
        UInteger rCount = context->argument(0).toUInt32();
        UInteger lCount = context->argument(1).toUInt32();
        double radius = context->argument(2).toNumber();
        double length = context->argument(3).toNumber();
        if (context->argumentCount() == 4)
        {
            QTriangleMesh3D *tmo = new QTriangleMesh3D();
            tmo->cylinderDomain(rCount, lCount, radius, length);

            return engine->newQObject(tmo, QScriptEngine::ScriptOwnership);
        }
        else
        {
            QScriptValue function = context->argument(4);
            if (!function.isFunction())
                return context->throwError(typeError.arg("func"));
            auto func = [&](double x, double y, double z)
            {
                QScriptValueList args;
                args << x << y << z;
                return function.call(QScriptValue(), args).toNumber();
            };

            QTriangleMesh3D *tmo = new QTriangleMesh3D();
            tmo->cylinderDomain(rCount, lCount, radius, length, func);

            return engine->newQObject(tmo, QScriptEngine::ScriptOwnership);
        }
    }
    return context->throwError(QObject::tr("CylinderTriangles(): arguments count error."));
}

QScriptValue QZScriptEngine::createConeTriangles(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() == 5 || context->argumentCount() == 6)
    {
        QString typeError = QObject::tr("ConeTriangles(rCount: Integer, lCount: Integer, radiusBottom: Floating, radiusTop: Floating, length: Floating): argument type error (%1).");
        if (!context->argument(0).isNumber())
            return context->throwError(typeError.arg("rCount"));
        if (!context->argument(1).isNumber())
            return context->throwError(typeError.arg("lCount"));
        if (!context->argument(2).isNumber())
            return context->throwError(typeError.arg("radiusBottom"));
        if (!context->argument(3).isNumber())
            return context->throwError(typeError.arg("radiusTop"));
        if (!context->argument(4).isNumber())
            return context->throwError(typeError.arg("length"));
        UInteger rCount = context->argument(0).toUInt32();
        UInteger lCount = context->argument(1).toUInt32();
        double radiusBottom = context->argument(2).toNumber();
        double radiusTop = context->argument(3).toNumber();
        double length = context->argument(4).toNumber();
        if (context->argumentCount() == 5)
        {
            QTriangleMesh3D *tmo = new QTriangleMesh3D();
            tmo->coneDomain(rCount, lCount, radiusBottom, radiusTop, length);
            return engine->newQObject(tmo, QScriptEngine::ScriptOwnership);
        }
        else
        {
            QScriptValue function = context->argument(5);
            if (!function.isFunction())
                return context->throwError(typeError.arg("func"));
            auto func = [&](double x, double y, double z)
            {
                QScriptValueList args;
                args << x << y << z;
                return function.call(QScriptValue(), args).toNumber();
            };
            QTriangleMesh3D *tmo = new QTriangleMesh3D();
            tmo->coneDomain(rCount, lCount, radiusBottom, radiusTop, length, func);
            return engine->newQObject(tmo, QScriptEngine::ScriptOwnership);
        }

    }
    return context->throwError(QObject::tr("ConeTriangles(): arguments count error."));
}

QScriptValue QZScriptEngine::createParametricTriangles(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() == 3 || context->argumentCount() == 4)
    {
        QString typeError = QObject::tr("ParametricTriangles(uCount: Integer, vCount: Integer, domain: Function [, rfunc: Function]): argument type error (%1).");
        if (!context->argument(0).isNumber())
            return context->throwError(typeError.arg("uCount"));
        if (!context->argument(1).isNumber())
            return context->throwError(typeError.arg("vCount"));
        QScriptValue domain_function = context->argument(2);
        if (!domain_function.isFunction())
            return context->throwError(typeError.arg("domain"));
        auto domain = [&](double u, double v)
        {
            QScriptValueList args;
            args << u << v;
            QPoint3D *p = qscriptvalue_cast<QPoint3D *>(domain_function.call(QScriptValue(), args));
            return Point3D(*p);
        };

        UInteger uCount = context->argument(0).toUInt32();
        UInteger vCount = context->argument(1).toUInt32();

        if (context->argumentCount() == 3)
        {
            QTriangleMesh3D *tmo = new QTriangleMesh3D();
            tmo->parametricDomain(uCount, vCount, domain, nullptr);
            return engine->newQObject(tmo, QScriptEngine::ScriptOwnership);
        }
        else
        {
            QScriptValue r_function = context->argument(3);
            if (!r_function.isFunction())
                return context->throwError(typeError.arg("rfunc"));
            auto rfunc = [&](double x, double y, double z)
            {
                QScriptValueList args;
                args << x << y << z;
                return (r_function.call(QScriptValue(), args)).toNumber();
            };
            QTriangleMesh3D *tmo = new QTriangleMesh3D();
            tmo->parametricDomain(uCount, vCount, domain, rfunc);
            return engine->newQObject(tmo, QScriptEngine::ScriptOwnership);
        }

    }
    return context->throwError(QObject::tr("ParametricTriangles(): arguments count error."));
}

QScriptValue QZScriptEngine::createMarchingCubes(QScriptContext *context, QScriptEngine *engine)
{
    if (8 <= context->argumentCount() && context->argumentCount() <= 11)
        {
            QString typeError = QObject::tr("MarchingCubes(xCount: Integer, yCount: Integer, zCount:Integer, origin: Point, width: Floating, height: Floating, depth: Floating, function: Function[, slice_x: Boolean, slice_y: Boolean, slice_z: Boolean]): argument type error (%1).");
            std::list<msh::Point2D> pointList;
            if (!context->argument(0).isNumber())
                return context->throwError(typeError.arg("xCount"));
            if (!context->argument(1).isNumber())
                return context->throwError(typeError.arg("yCount"));
            if (!context->argument(2).isNumber())
                return context->throwError(typeError.arg("zCount"));
            if (!context->argument(3).isQObject() || qscriptvalue_cast<QPoint3D *>(context->argument(3)) == NULL)
                return context->throwError(typeError.arg("origin"));
            if (!context->argument(4).isNumber())
                return context->throwError(typeError.arg("width"));
            if (!context->argument(5).isNumber())
                return context->throwError(typeError.arg("height"));
            if (!context->argument(6).isNumber())
                return context->throwError(typeError.arg("depth"));
            QScriptValue function = context->argument(7);
            if (!function.isFunction())
                return context->throwError(typeError.arg("function"));
            if (context->argumentCount() >= 9 && !context->argument(8).isBool())
                return context->throwError(typeError.arg("slice_x"));
            if (context->argumentCount() >= 10 && !context->argument(9).isBool())
                return context->throwError(typeError.arg("slice_y"));
            if (context->argumentCount() == 11 && !context->argument(10).isBool())
                return context->throwError(typeError.arg("slice_z"));
            UInteger xCount = context->argument(0).toUInt32();
            UInteger yCount = context->argument(1).toUInt32();
            UInteger zCount = context->argument(2).toUInt32();
            QPoint3D *origin = qscriptvalue_cast<QPoint3D *>(context->argument(3));
            double width = context->argument(4).toNumber();
            double height = context->argument(5).toNumber();
            double depth = context->argument(6).toNumber();
            bool slice_x = (context->argumentCount() >= 9) ? context->argument(8).toBool() : false;
            bool slice_y = (context->argumentCount() >= 10) ? context->argument(9).toBool() : false;
            bool slice_z = (context->argumentCount() == 11) ? context->argument(10).toBool() : false;
            // R-функция
            auto func = [&](double x, double y, double z)
            {
                QScriptValueList args;
                args << x << y << z;
                return function.call(QScriptValue(), args).toNumber();
            };

            QTriangleMesh3D *tmo = new QTriangleMesh3D();
            tmo->marchingCubes(xCount, yCount, zCount, origin->x(), origin->y(), origin->z(), width, height, depth, func, 0.0, slice_x, slice_y, slice_z);

            return engine->newQObject(tmo, QScriptEngine::ScriptOwnership);
        }
        return context->throwError(QObject::tr("MarchingCubes(xCount: Integer, yCount: Integer, zCount:Integer, origin: Point, width: Floating, height: Floating, depth: Floating, function: Function): arguments count error."));
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
    else if (qscriptvalue_cast<QSegmentMesh2D *>(context->argument(0)) != NULL)
    {
        mesh_ = new SegmentMesh2D(qscriptvalue_cast<QSegmentMesh2D *>(context->argument(0)));
    }
    else if (qscriptvalue_cast<QTriangleMesh2D *>(context->argument(0)) != NULL)
    {
        mesh_ = new TriangleMesh2D(qscriptvalue_cast<QTriangleMesh2D *>(context->argument(0)));
    }
    else if (qscriptvalue_cast<QQuadrilateralMesh3D *>(context->argument(0)) != NULL)
    {
        mesh_ = new QuadrilateralMesh3D(qscriptvalue_cast<QQuadrilateralMesh3D *>(context->argument(0)));
    }
    else if (qscriptvalue_cast<QTriangleMesh3D *>(context->argument(0)) != NULL)
    {
        mesh_ = new TriangleMesh3D(qscriptvalue_cast<QTriangleMesh3D *>(context->argument(0)));
    }
    else
    {
        return context->throwError(typeError);
    }
    return engine->undefinedValue();
}

QScriptValue QZScriptEngine::currentMesh(QScriptContext *context, QScriptEngine *engine)
{
    if (dynamic_cast<msh::SegmentMesh2D *>(mesh_))
        return engine->newQObject(new QSegmentMesh2D(dynamic_cast<SegmentMesh2D *>(mesh_)), QScriptEngine::ScriptOwnership);
    else if (dynamic_cast<msh::TriangleMesh2D *>(mesh_))
        return engine->newQObject(new QTriangleMesh2D(dynamic_cast<msh::TriangleMesh2D *>(mesh_)), QScriptEngine::ScriptOwnership);
    else if (dynamic_cast<msh::QuadrilateralMesh2D *>(mesh_))
        return engine->newQObject(new QQuadrilateralMesh2D(dynamic_cast<msh::QuadrilateralMesh2D *>(mesh_)), QScriptEngine::ScriptOwnership);
    else if (dynamic_cast<msh::TriangleMesh3D *>(mesh_))
        return engine->newQObject(new QTriangleMesh3D(dynamic_cast<msh::TriangleMesh3D *>(mesh_)), QScriptEngine::ScriptOwnership);
    else if (dynamic_cast<msh::QuadrilateralMesh3D *>(mesh_))
        return engine->newQObject(new QQuadrilateralMesh3D(dynamic_cast<msh::QuadrilateralMesh3D *>(mesh_)), QScriptEngine::ScriptOwnership);
    else
        return context->throwError(QObject::tr("Type of current mesh is unusefull."));
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
                    return context->throwError(typeError.arg("Point"));
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

QScriptValue QZScriptEngine::circle(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);
    QString typeError = QObject::tr("circle(x, y, r): all arguments must have type Number. Argument # %1.");
    if (context->argumentCount() != 3)
        return context->throwError(tr("Function circle has at least three arguments."));
    for (int i = 0; i < context->argumentCount(); i++)
        if (!context->argument(i).isNumber())
            return context->throwError(typeError.arg(i + 1));
    double x = context->argument(0).toNumber();
    double y = context->argument(1).toNumber();
    double r = context->argument(2).toNumber();
    return msh::circle(x, y, r);
}

QScriptValue QZScriptEngine::ellipse(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);
    QString typeError = QObject::tr("ellipse(x, y, a, b): all arguments must have type Number. Argument # %1.");
    if (context->argumentCount() != 4)
        return context->throwError(tr("Function circle has at least four arguments."));
    for (int i = 0; i < context->argumentCount(); i++)
        if (!context->argument(i).isNumber())
            return context->throwError(typeError.arg(i + 1));
    double x = context->argument(0).toNumber();
    double y = context->argument(1).toNumber();
    double a = context->argument(2).toNumber();
    double b = context->argument(3).toNumber();
    return msh::ellipse(x, y, a, b);
}

QScriptValue QZScriptEngine::band(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);
    QString typeError = QObject::tr("band(x, w): all arguments must have type Number. Argument # %1.");
    if (context->argumentCount() != 2)
        return context->throwError(tr("Function circle has at least two arguments."));
    for (int i = 0; i < context->argumentCount(); i++)
        if (!context->argument(i).isNumber())
            return context->throwError(typeError.arg(i + 1));
    double x = context->argument(0).toNumber();
    double w = context->argument(1).toNumber();
    return msh::band(x, w);
}

QScriptValue QZScriptEngine::line(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);
    QString typeError = QObject::tr("band(x, w): all arguments must have type Number. Argument # %1.");
    if (context->argumentCount() != 6)
        return context->throwError(tr("Function circle has at least six arguments."));
    for (int i = 0; i < context->argumentCount(); i++)
        if (!context->argument(i).isNumber())
            return context->throwError(typeError.arg(i + 1));
    double x = context->argument(0).toNumber();
    double y = context->argument(1).toNumber();
    double x1 = context->argument(2).toNumber();
    double y1 = context->argument(3).toNumber();
    double x2 = context->argument(4).toNumber();
    double y2 = context->argument(5).toNumber();
    return msh::line(x, y, x1, y1, x2, y2);
}

QScriptValue QZScriptEngine::rectangle(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);
    QString typeError = QObject::tr("rectangle(x, y, w, h [, r]): all arguments must have type Number. Argument # %1.");
    if (context->argumentCount() != 4 && context->argumentCount() != 5)
        return context->throwError(tr("Function rectangle has at least four arguments."));
    for (int i = 0; i < context->argumentCount(); i++)
        if (!context->argument(i).isNumber())
            return context->throwError(typeError.arg(i + 1));
    double x = context->argument(0).toNumber();
    double y = context->argument(1).toNumber();
    double w = context->argument(2).toNumber();
    double h = context->argument(3).toNumber();
    if (context->argumentCount() == 4)
        return msh::rectangle(x, y, w, h);
    double r = context->argument(4).toNumber();
    if (r < 0.0)
        return context->throwError(typeError.arg("r < 0"));
    return msh::rectangle(x, y, w, h, r);
}

QScriptValue QZScriptEngine::convex(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);
    QString typeError = QObject::tr("convex(x, y, P: Points): argument type error (%1).");
    if (context->argumentCount() < 3)
        context->throwError(QObject::tr("convex(x, y, P: Points): arguments count error."));
    if (!context->argument(0).isNumber())
        return context->throwError(typeError.arg("x: Float"));
    if (!context->argument(1).isNumber())
        return context->throwError(typeError.arg("y: Float"));
    double x = context->argument(0).toNumber();
    double y = context->argument(1).toNumber();
    if (context->argument(2).isArray())
    {
        QScriptValue array = context->argument(2);
        double r = 0.0;
        int pcount = array.property("length").toInteger();
        if (pcount < 3)
            return context->throwError(typeError.arg("P is array of at least 3 points."));
        for (int i = 0; i < pcount; i++)
        {
            if (qscriptvalue_cast<QPoint2D *>(array.property(i)) == NULL)
                return context->throwError(typeError.arg(QString("Broken point value # ") + QString::number(i)));
        }
        QPoint2D *p1 = qscriptvalue_cast<QPoint2D*>(array.property(0));
        QPoint2D *p2 = qscriptvalue_cast<QPoint2D*>(array.property(1));
        r = msh::line(x, y, p1->x(), p1->y(), p2->x(), p2->y());
        for (int i = 2; i < pcount; i++)
        {
            p1 = p2;
            p2 = qscriptvalue_cast<QPoint2D*>(array.property(i));
            r = msh::con(r, msh::line(x, y, p1->x(), p1->y(), p2->x(), p2->y()));
        }
        p1 = qscriptvalue_cast<QPoint2D*>(array.property(pcount - 1));
        p2 = qscriptvalue_cast<QPoint2D*>(array.property(0));
        r = msh::con(r, msh::line(x, y, p1->x(), p1->y(), p2->x(), p2->y()));
        return r;
    }
    return context->throwError(typeError.arg("p: Points"));

}

QScriptValue QZScriptEngine::regular(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);
    QString typeError = QObject::tr("regular(x, y, r: Float, n: Integer): argument type error (%1).");
    if (context->argumentCount() != 4)
        context->throwError(QObject::tr("regular(x, y, r: Float, n: Integer): arguments count error."));
    if (!context->argument(0).isNumber())
        return context->throwError(typeError.arg("x: Floating number"));
    if (!context->argument(1).isNumber())
        return context->throwError(typeError.arg("y: Floating number"));
    if (!context->argument(2).isNumber())
        return context->throwError(typeError.arg("r: Floating number"));
    if (!context->argument(3).isNumber())
        return context->throwError(typeError.arg("n: Integer number"));
    double x = context->argument(0).toNumber();
    double y = context->argument(1).toNumber();
    double r = context->argument(2).toNumber();
    int n = context->argument(3).toInt32();
    return msh::regular(x, y, r, n);
}

QScriptValue QZScriptEngine::ellipsoid(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);
    QString typeError = QObject::tr("ellipsoid(x, y, z, a, b, c: Float): argument type error (%1).");
    if (context->argumentCount() != 6)
        context->throwError(QObject::tr("ellipsoid(x, y, z, a, b, c): arguments count error."));
    if (!context->argument(0).isNumber())
        return context->throwError(typeError.arg("x: Floating number"));
    if (!context->argument(1).isNumber())
        return context->throwError(typeError.arg("y: Floating number"));
    if (!context->argument(2).isNumber())
        return context->throwError(typeError.arg("z: Floating number"));
    if (!context->argument(3).isNumber())
        return context->throwError(typeError.arg("a: Floating number"));
    if (!context->argument(4).isNumber())
        return context->throwError(typeError.arg("b: Floating number"));
    if (!context->argument(5).isNumber())
        return context->throwError(typeError.arg("c: Floating number"));
    double x = context->argument(0).toNumber();
    double y = context->argument(1).toNumber();
    double z = context->argument(2).toNumber();
    double a = context->argument(3).toNumber();
    double b = context->argument(4).toNumber();
    double c = context->argument(5).toNumber();
    return msh::ellipsoid(x, y, z, a, b, c);
}

QScriptValue QZScriptEngine::sphere(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);
    QString typeError = QObject::tr("sphere(x, y, z, r: Float): argument type error (%1).");
    if (context->argumentCount() != 4)
        context->throwError(QObject::tr("sphere(x, y, z, r): arguments count error."));
    if (!context->argument(0).isNumber())
        return context->throwError(typeError.arg("x: Floating number"));
    if (!context->argument(1).isNumber())
        return context->throwError(typeError.arg("y: Floating number"));
    if (!context->argument(2).isNumber())
        return context->throwError(typeError.arg("z: Floating number"));
    if (!context->argument(3).isNumber())
        return context->throwError(typeError.arg("r: Floating number"));
    double x = context->argument(0).toNumber();
    double y = context->argument(1).toNumber();
    double z = context->argument(2).toNumber();
    double r = context->argument(3).toNumber();
    return msh::sphere(x, y, z, r);
}

QScriptValue QZScriptEngine::plane(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);
    QString typeError = QObject::tr("plane(x, y, z: Float, P1, P2, P3: Point): argument type error (%1).");
    if (context->argumentCount() != 6)
        context->throwError(QObject::tr("plane(x, y, z, P1, P2, P3): arguments count error."));
    if (!context->argument(0).isNumber())
        return context->throwError(typeError.arg("x: Floating number"));
    if (!context->argument(1).isNumber())
        return context->throwError(typeError.arg("y: Floating number"));
    if (!context->argument(2).isNumber())
        return context->throwError(typeError.arg("z: Floating number"));
    QPoint3D *p1 = qscriptvalue_cast<QPoint3D*>(context->argument(3));
    QPoint3D *p2 = qscriptvalue_cast<QPoint3D*>(context->argument(4));
    QPoint3D *p3 = qscriptvalue_cast<QPoint3D*>(context->argument(5));
    if (p1 == NULL)
        return context->throwError(typeError.arg("p1: 3d point (Point)"));
    if (p2 == NULL)
        return context->throwError(typeError.arg("p2: 3d point (Point)"));
    if (p3 == NULL)
        return context->throwError(typeError.arg("p3: 3d point (Point)"));
    double x = context->argument(0).toNumber();
    double y = context->argument(1).toNumber();
    double z = context->argument(2).toNumber();
    return msh::plane(x, y, z, p1->x(), p1->y(), p1->z(), p2->x(), p2->y(), p2->z(), p3->x(), p3->y(), p3->z());
}

QScriptValue QZScriptEngine::cuboid(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);
    QString typeError = QObject::tr("cuboid(x, y, z, width, height, depth: Float): argument type error (%1).");
    if (context->argumentCount() != 6)
        context->throwError(QObject::tr("cuboid(x, y, z, width, height, depth): arguments count error."));
    if (!context->argument(0).isNumber())
        return context->throwError(typeError.arg("x: Floating number"));
    if (!context->argument(1).isNumber())
        return context->throwError(typeError.arg("y: Floating number"));
    if (!context->argument(2).isNumber())
        return context->throwError(typeError.arg("z: Floating number"));
    if (!context->argument(3).isNumber())
        return context->throwError(typeError.arg("width: Floating number"));
    if (!context->argument(4).isNumber())
        return context->throwError(typeError.arg("height: Floating number"));
    if (!context->argument(5).isNumber())
        return context->throwError(typeError.arg("delpth: Floating number"));
    double x = context->argument(0).toNumber();
    double y = context->argument(1).toNumber();
    double z = context->argument(2).toNumber();
    double w = context->argument(3).toNumber();
    double h = context->argument(4).toNumber();
    double d = context->argument(5).toNumber();
    return msh::cuboid(x, y, z, w, h, d);
}

QScriptValue QZScriptEngine::planeStress(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() > 4)
    {
        Mesh2D *mesh = NULL;
        QString typeError = QObject::tr("PlaneStress(mesh: Mesh, h: Floating, E: Floating, nu: Floating, ...): argument type error (%1).");
        if (!context->argument(0).isQObject())
            return context->throwError(typeError.arg("mesh"));
        if (qscriptvalue_cast<QQuadrilateralMesh2D *>(context->argument(0)) != NULL)
        {
            mesh = new QuadrilateralMesh2D(qscriptvalue_cast<QQuadrilateralMesh2D *>(context->argument(0)));
        }
        else if (qscriptvalue_cast<QTriangleMesh2D *>(context->argument(0)) != NULL)
        {
            mesh = new TriangleMesh2D(qscriptvalue_cast<QTriangleMesh2D *>(context->argument(0)));
        }
        else
        {
            return context->throwError(typeError.arg("mesh"));
        }

        if (!context->argument(1).isNumber())
            return context->throwError(typeError.arg("h"));

        if (!context->argument(2).isNumber())
            return context->throwError(typeError.arg("E"));

        if (!context->argument(3).isNumber())
            return context->throwError(typeError.arg("nu"));

        std::list<FemCondition*> conditions;

        for (int i = 4; i < context->argumentCount(); i++)
        {
            if (!context->argument(i).isQObject())
                return context->throwError(typeError.arg("boundary condition"));
            QFemCondition *cond = qscriptvalue_cast<QFemCondition *>(context->argument(i));
            if (cond == NULL)
                return context->throwError(typeError.arg("boundary condition"));
            conditions.push_back(cond);
        }

        double h = context->argument(1).toNumber();
        double E = context->argument(2).toNumber();
        double nu = context->argument(3).toNumber();

        DoubleMatrix D = Fem2D::evalPlaneStressMatrix(E, nu);

        if (fem_ != NULL) delete fem_;

        fem_ = new PlaneStressStrain (mesh, //!
                                      h,
                                      D,
                                      conditions);
        fem_->solve();
        mesh->printDataExtremums();

        setMesh(mesh);

        return engine->undefinedValue();;
    }
    return context->throwError(QObject::tr("PlaneStress(mesh: Mesh, h: Floating, E: Floating, nu: Floating, boundaryConditionType: Function, boundaryValue: Function, nodalForce: Function, surfaceForce: Function, volumeForce: Function): arguments count error."));
}

QScriptValue QZScriptEngine::planeStrain(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() > 4)
    {
        Mesh2D *mesh = NULL;
        QString typeError = QObject::tr("PlaneStrain(mesh: Mesh, h: Floating, E: Floating, nu: Floating, ...): argument type error (%1).");
        if (!context->argument(0).isQObject())
            return context->throwError(typeError.arg("mesh"));
        if (qscriptvalue_cast<QQuadrilateralMesh2D *>(context->argument(0)) != NULL)
        {
            mesh = new QuadrilateralMesh2D(qscriptvalue_cast<QQuadrilateralMesh2D *>(context->argument(0)));
        }
        else if (qscriptvalue_cast<QTriangleMesh2D *>(context->argument(0)) != NULL)
        {
            mesh = new TriangleMesh2D(qscriptvalue_cast<QTriangleMesh2D *>(context->argument(0)));
        }
        else
        {
            return context->throwError(typeError.arg("mesh"));
        }

        if (!context->argument(1).isNumber())
            return context->throwError(typeError.arg("h"));

        if (!context->argument(2).isNumber())
            return context->throwError(typeError.arg("E"));

        if (!context->argument(3).isNumber())
            return context->throwError(typeError.arg("nu"));

        std::list<FemCondition*> conditions;

        for (int i = 4; i < context->argumentCount(); i++)
        {
            if (!context->argument(i).isQObject())
                return context->throwError(typeError.arg("boundary condition"));
            QFemCondition *cond = qscriptvalue_cast<QFemCondition *>(context->argument(i));
            if (cond == NULL)
                return context->throwError(typeError.arg("boundary condition"));
            conditions.push_back(cond);
        }

        double h = context->argument(1).toNumber();
        double E = context->argument(2).toNumber();
        double nu = context->argument(3).toNumber();

        DoubleMatrix D = Fem2D::evalPlaneStrainMatrix(E, nu); //!

        if (fem_ != NULL) delete fem_;

        fem_ = new PlaneStressStrain (mesh, //!
                                      h,
                                      D,
                                      conditions);
        fem_->solve();
        mesh->printDataExtremums();

        setMesh(mesh);

        return engine->undefinedValue();;
    }
    return context->throwError(QObject::tr("PlaneStrain(mesh: Mesh, h: Floating, E: Floating, nu: Floating, boundaryConditionType: Function, boundaryValue: Function, nodalForce: Function, surfaceForce: Function, volumeForce: Function): arguments count error."));
}

QScriptValue QZScriptEngine::createBoundaryCondition(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() == 3)
    {
        QString typeError = QObject::tr("BoundaryCondition(direction: Integer, condition: Function, value: {Floating or Function}): argument type error (%1)");
        if (!context->argument(0).isNumber())
            return context->throwError(typeError.arg("direction"));

        QScriptValue condition = context->argument(1);
        if (!condition.isFunction())
            return context->throwError(typeError.arg("condition"));

        QScriptValue value = context->argument(2);

        int direction = context->argument(0).toInteger();

        return engine->newQObject(new QFemCondition(FemCondition::INITIAL_VALUE, static_cast<FemCondition::FemDirection>(direction), condition, value), QScriptEngine::ScriptOwnership);
    }
    return context->throwError(QObject::tr("BoundaryCondition(direction: Integer, condition: Function, value: {Floating or Function}): arguments count error."));
}

QScriptValue QZScriptEngine::createNodalForce(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() == 2)
    {   // Если аргумента 2, то сила действует на все точки
        QString typeError = QObject::tr("NodalForce(direction: Integer, value: {Floating or Function}): argument type error (%1)");
        if (!context->argument(0).isNumber())
            return context->throwError(typeError.arg("direction"));

        QScriptValue value = context->argument(1);

        int direction = context->argument(0).toInteger();

        return engine->newQObject(new QFemCondition(FemCondition::NODAL_FORCE, static_cast<FemCondition::FemDirection>(direction), QScriptValue(true), value), QScriptEngine::ScriptOwnership);
    }
    if (context->argumentCount() == 3)
    {
        QString typeError = QObject::tr("NodalForce(direction: Integer, condition: Function, value: {Floating or Function}): argument type error (%1)");
        if (!context->argument(0).isNumber())
            return context->throwError(typeError.arg("direction"));

        QScriptValue condition = context->argument(1);
        if (!condition.isFunction())
            return context->throwError(typeError.arg("condition"));

        QScriptValue value = context->argument(2);

        int direction = context->argument(0).toInteger();

        return engine->newQObject(new QFemCondition(FemCondition::NODAL_FORCE, static_cast<FemCondition::FemDirection>(direction), condition, value), QScriptEngine::ScriptOwnership);
    }
    return context->throwError(QObject::tr("NodalForce(direction: Integer, condition: Function, value: {Floating or Function}): arguments count error."));
}

QScriptValue QZScriptEngine::createSurfaceForce(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() == 2)
    {   // Если аргумента 2, то сила действует на все точки
        QString typeError = QObject::tr("SurfaceForce(direction: Integer, value: {Floating or Function}): argument type error (%1)");
        if (!context->argument(0).isNumber())
            return context->throwError(typeError.arg("direction"));

        QScriptValue value = context->argument(1);

        int direction = context->argument(0).toInteger();

        return engine->newQObject(new QFemCondition(FemCondition::SURFACE_FORCE, static_cast<FemCondition::FemDirection>(direction), QScriptValue(true), value), QScriptEngine::ScriptOwnership);
    }
    if (context->argumentCount() == 3)
    {
        QString typeError = QObject::tr("SurfaceForce(direction: Integer, condition: Function, value: {Floating or Function}): argument type error (%1)");
        if (!context->argument(0).isNumber())
            return context->throwError(typeError.arg("direction"));

        QScriptValue condition = context->argument(1);
        if (!condition.isFunction())
            return context->throwError(typeError.arg("condition"));

        QScriptValue value = context->argument(2);

        int direction = context->argument(0).toInteger();

        return engine->newQObject(new QFemCondition(FemCondition::SURFACE_FORCE, static_cast<FemCondition::FemDirection>(direction), condition, value), QScriptEngine::ScriptOwnership);
    }
    return context->throwError(QObject::tr("SurfaceForce(direction: Integer, condition: Function, value: {Floating or Function}): arguments count error."));
}

QScriptValue QZScriptEngine::createVolumeForce(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() == 2)
    {   // Если аргумента 2, то сила действует на все точки
        QString typeError = QObject::tr("VolumeForce(direction: Integer, value: {Floating or Function}): argument type error (%1)");
        if (!context->argument(0).isNumber())
            return context->throwError(typeError.arg("direction"));

        QScriptValue value = context->argument(1);

        int direction = context->argument(0).toInteger();

        return engine->newQObject(new QFemCondition(FemCondition::VOLUME_FORCE, static_cast<FemCondition::FemDirection>(direction), QScriptValue(true), value), QScriptEngine::ScriptOwnership);
    }
    if (context->argumentCount() == 3)
    {
        QString typeError = QObject::tr("VolumeForce(direction: Integer, condition: Function, value: {Floating or Function}): argument type error (%1)");
        if (!context->argument(0).isNumber())
            return context->throwError(typeError.arg("direction"));

        QScriptValue condition = context->argument(1);
        if (!condition.isFunction())
            return context->throwError(typeError.arg("condition"));

        QScriptValue value = context->argument(2);

        int direction = context->argument(0).toInteger();

        return engine->newQObject(new QFemCondition(FemCondition::VOLUME_FORCE, static_cast<FemCondition::FemDirection>(direction), condition, value), QScriptEngine::ScriptOwnership);
    }
    return context->throwError(QObject::tr("VolumeForce(direction: Integer, condition: Function, value: {Floating or Function}): arguments count error."));
}

QScriptValue QZScriptEngine::mindlinPlate(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() > 4)
    {
        Mesh2D *mesh = NULL;
        QString typeError = QObject::tr("MindlinPlate(mesh: Mesh, h: Floating, E: Floating, nu: Floating, ...): argument type error (%1).");
        if (!context->argument(0).isQObject())
            return context->throwError(typeError.arg("mesh"));
        if (qscriptvalue_cast<QQuadrilateralMesh2D *>(context->argument(0)) != NULL)
        {
            mesh = new QuadrilateralMesh2D(qscriptvalue_cast<QQuadrilateralMesh2D *>(context->argument(0)));
        }
        else if (qscriptvalue_cast<QTriangleMesh2D *>(context->argument(0)) != NULL)
        {
            mesh = new TriangleMesh2D(qscriptvalue_cast<QTriangleMesh2D *>(context->argument(0)));
        }
        else
        {
            return context->throwError(typeError.arg("mesh"));
        }

        if (!context->argument(1).isNumber() && !context->argument(1).isArray())
            return context->throwError(typeError.arg("h"));

        if (!context->argument(2).isNumber() && !context->argument(2).isArray())
            return context->throwError(typeError.arg("E"));

        if (!context->argument(3).isNumber() && !context->argument(3).isArray())
            return context->throwError(typeError.arg("nu"));

        std::list<FemCondition*> conditions;

        for (int i = 4; i < context->argumentCount(); i++)
        {
            if (!context->argument(i).isQObject())
                return context->throwError(typeError.arg("boundary condition"));
            QFemCondition *cond = qscriptvalue_cast<QFemCondition *>(context->argument(i));
            if (cond == NULL)
                return context->throwError(typeError.arg("boundary condition"));
            conditions.push_back(cond);
        }
        if (context->argument(1).isNumber() && context->argument(2).isNumber() && context->argument(3).isNumber())
        {
            double h = context->argument(1).toNumber();
            double E = context->argument(2).toNumber();
            double nu = context->argument(3).toNumber();

            DoubleMatrix D = Fem2D::evalPlaneStressMatrix(E, nu); //!

            if (fem_ != NULL) delete fem_;

            fem_ = new MindlinPlateBending (mesh, //!
                                            h,
                                            D,
                                            conditions);
            fem_->solve();
            setMesh(mesh);
        }
        else if (context->argument(1).isArray() && context->argument(2).isArray() && context->argument(3).isArray())
        {
            QScriptValue h_array = context->argument(1);
            QScriptValue e_array = context->argument(2);
            QScriptValue nu_array = context->argument(3);
            std::vector<double> h;
            std::vector<DoubleMatrix> elasticMatrix;

            if (h_array.property("length").toInteger() != e_array.property("length").toInteger() || h_array.property("length").toInteger() != nu_array.property("length").toInteger())
            {
                return context->throwError(typeError.arg("h, E, nu: all the input arrays must have same length"));
            }

            for (int i = 0; i < h_array.property("length").toInteger(); i++)
            {
                h.push_back(h_array.property(i).toNumber());
                elasticMatrix.push_back(Fem2D::evalPlaneStressMatrix(e_array.property(i).toNumber(), nu_array.property(i).toNumber()));
            }

            if (fem_ != NULL) delete fem_;

            fem_ = new MindlinPlateLaminated (mesh, //!
                                            h,
                                            elasticMatrix,
                                            conditions);
            fem_->solve();
            setMesh(mesh);
        }
        else
        {
            return context->throwError(typeError.arg("h, E, nu: all the input must be scalars or arrays of same length"));
        }

        mesh->printDataExtremums();

        return engine->undefinedValue();
    }
    return context->throwError(QObject::tr("MindlinPlate(mesh: Mesh, h: Floating, E: Floating, nu: Floating, boundaryConditionType: Function, boundaryValue: Function, nodalForce: Function, surfaceForce: Function, volumeForce: Function): arguments count error."));
}

QScriptValue QZScriptEngine::mindlinShell(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() > 4)
    {
        Mesh3D *mesh = NULL;
        QString typeError = QObject::tr("MindlinShell(mesh: Mesh, h: Floating, {E: Floating, nu: Floating[, G: Floating]} | {strain: Array, stress: Array, nu: Floating}, [, {boundary conditions, forces}]...): argument type error (%1).");
        if (!context->argument(0).isQObject())
            return context->throwError(typeError.arg("mesh"));
        if (qscriptvalue_cast<QQuadrilateralMesh3D *>(context->argument(0)) != NULL)
        {
            mesh = new QuadrilateralMesh3D(qscriptvalue_cast<QQuadrilateralMesh3D *>(context->argument(0)));
        }
        else if (qscriptvalue_cast<QTriangleMesh3D *>(context->argument(0)) != NULL)
        {
            mesh = new TriangleMesh3D(qscriptvalue_cast<QTriangleMesh3D *>(context->argument(0)));
        }
        else
        {
            return context->throwError(typeError.arg("mesh"));
        }

        if (!context->argument(1).isNumber() && !context->argument(1).isArray())
            return context->throwError(typeError.arg("h"));

        if (!context->argument(2).isNumber() && !context->argument(2).isArray())
            return context->throwError(typeError.arg("E"));

        if (!context->argument(3).isNumber() && !context->argument(3).isArray())
            return context->throwError(typeError.arg("nu"));

        int start_pos = 4;
        if (!context->argument(4).isNumber() && !context->argument(4).isArray())
            start_pos = 4;
        else
            start_pos = 5;

        std::list<FemCondition*> conditions;

        for (int i = start_pos; i < context->argumentCount(); i++)
        {
            if (!context->argument(i).isQObject())
                return context->throwError(typeError.arg("boundary condition or force"));
            QFemCondition *cond = qscriptvalue_cast<QFemCondition *>(context->argument(i));
            if (cond == NULL)
                return context->throwError(typeError.arg("boundary condition or force"));
            conditions.push_back(cond);
        }
        if (context->argument(1).isNumber() && context->argument(2).isNumber() && context->argument(3).isNumber())
        {
            double h = context->argument(1).toNumber();
            double E = context->argument(2).toNumber();
            double nu = context->argument(3).toNumber();

            DoubleMatrix D;

            if (context->argument(4).isNumber())
                D = Fem2D::evalPlaneStressMatrix(E, nu, context->argument(4).toNumber());
            else
                D = Fem2D::evalPlaneStressMatrix(E, nu);

            if (fem_ != NULL) delete fem_;

            fem_ = new MindlinShellBending (mesh, //!
                                            h,
                                            D,
                                            conditions);
            fem_->solve();
            setMesh(mesh);
        }
        else if (context->argument(1).isArray() && context->argument(2).isArray() && context->argument(3).isArray())
        {
            QScriptValue h_array = context->argument(1);
            QScriptValue e_array = context->argument(2);
            QScriptValue nu_array = context->argument(3);
            std::vector<double> h;
            std::vector<DoubleMatrix> elasticMatrix;

            if (h_array.property("length").toInteger() != e_array.property("length").toInteger() || h_array.property("length").toInteger() != nu_array.property("length").toInteger())
            {
                return context->throwError(typeError.arg("h, E, nu: all the input arrays must have same length"));
            }
            if (context->argument(4).isArray())
            {
                QScriptValue g_array = context->argument(4);
                if (h_array.property("length").toInteger() != g_array.property("length").toInteger())
                {
                    return context->throwError(typeError.arg("h, E, nu, G: all the input arrays must have same length"));
                }
                for (int i = 0; i < h_array.property("length").toInteger(); i++)
                {
                    DoubleMatrix D = Fem2D::evalPlaneStressMatrix(e_array.property(i).toNumber(), nu_array.property(i).toNumber(), g_array.property(i).toNumber());
                    h.push_back(h_array.property(i).toNumber());
                    elasticMatrix.push_back(D);
                }
            }
            else
            {
                for (int i = 0; i < h_array.property("length").toInteger(); i++)
                {
                    DoubleMatrix D = Fem2D::evalPlaneStressMatrix(e_array.property(i).toNumber(), nu_array.property(i).toNumber());
                    h.push_back(h_array.property(i).toNumber());
                    elasticMatrix.push_back(D);
                }
            }

            if (fem_ != NULL) delete fem_;

            fem_ = new MindlinShellLaminated (mesh, //!
                                            h,
                                            elasticMatrix,
                                            conditions);
            fem_->solve();
            setMesh(mesh);
        }
        else if (context->argument(1).isNumber() && context->argument(2).isArray() && context->argument(3).isArray() && context->argument(4).isNumber())
        {
            double h = context->argument(1).toNumber();
            QScriptValue strain_array = context->argument(2);
            QScriptValue stress_array = context->argument(3);
            double nu = context->argument(4).toNumber();
            std::vector<double> strain;
            std::vector<double> stress;
            if (strain_array.property("length").toInteger() != stress_array.property("length").toInteger())
            {
                return context->throwError(typeError.arg("stress, strain: all the input arrays must have same length"));
            }
            for (int i = 0; i < strain_array.property("length").toInteger(); i++)
            {
                strain.push_back(strain_array.property(i).toNumber());
                stress.push_back(stress_array.property(i).toNumber());
            }

            if (fem_ != NULL) delete fem_;

//            fem_ = new MindlinShellBending (mesh, h, strain, stress, nu, conditions);

            setMesh(mesh);
        }
        else
        {
            return context->throwError(typeError.arg("h, E, nu: all the input must be scalars or arrays of same length"));
        }

        mesh->printDataExtremums();

        return engine->undefinedValue();
    }
    return context->throwError(QObject::tr("MindlinShell(mesh: Mesh, h: Floating, E: Floating, nu: Floating, boundaryConditionType: Function, boundaryValue: Function, nodalForce: Function, surfaceForce: Function, volumeForce: Function): arguments count error."));
}

QScriptValue QZScriptEngine::reportValues(QScriptContext *context, QScriptEngine *engine)
{
     if (context->argumentCount() == 1 && fem_ != NULL)
     {
         QString typeError = QObject::tr("values(func: Function): argument type error (%1).");
         QScriptValue function = context->argument(0);
         if (!function.isFunction())
             return context->throwError(typeError.arg("func"));
         auto func = [&](msh::PointPointer point)
         {
             QScriptValueList args;
             args << point->x() << point->y() << point->z();
             return function.call(QScriptValue(), args).toBool();
         };
         fem_->reportNodeValues(func);
         return engine->undefinedValue();
     }
     return context->throwError(QObject::tr("values(func: Function): arguments count error."));
}

