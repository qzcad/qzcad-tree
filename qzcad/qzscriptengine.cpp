#include <iostream>
#include <math.h>
#include <list>
#include "qpoint2d.h"
#include "qpoint3d.h"
#include "qquadrilateralmesh2d.h"
#include "qtrianglemesh2d.h"
#include "qquadrilateralmesh3d.h"

#include "qfemcondition.h"

#include "planestressstrain.h"
#include "mindlinplatebending.h"

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
    // Поверхностная сетка четырехугольников
    QScriptValue qsCreateQuadrilateralMesh3D = newFunction(createQuadrilateralMesh3D);
    globalObject().setProperty("SurfaceQuads", qsCreateQuadrilateralMesh3D);;
    // setMesh
    QScriptValue qsSetMesh = newFunction(setMesh);
    globalObject().setProperty("setMesh", qsSetMesh);
    // About
    QScriptValue qsAbout = newFunction(about);
    globalObject().setProperty("About", qsAbout);
    /// FEM section
    QScriptValue qsPlaneStress = newFunction(planeStress);
    globalObject().setProperty("PlaneStress", qsPlaneStress);
    QScriptValue qsPlaneStrain = newFunction(planeStrain);
    globalObject().setProperty("PlaneStrain", qsPlaneStrain);

    QScriptValue qsBoundaryCondition = newFunction(createBoundaryCondition);
    globalObject().setProperty("BoundaryCondition", qsBoundaryCondition);
    QScriptValue qsNodalForce = newFunction(createNodalForce);
    globalObject().setProperty("NodalForce", qsNodalForce);
    QScriptValue qsSurfaceForce = newFunction(createSurfaceForce);
    globalObject().setProperty("SurfaceForce", qsSurfaceForce);
    QScriptValue qsVolumeForce = newFunction(createVolumeForce);
    globalObject().setProperty("VolumeForce", qsVolumeForce);
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

    QScriptValue qsMindlinPlate = newFunction(mindlinPlate);
    globalObject().setProperty("MindlinPlate", qsMindlinPlate);

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
    if (fem_ == NULL) return 0;
    return fem_->nodesVectorsCount();
}

NamedFloatingVector QZScriptEngine::getNodeValues(const unsigned long &i)
{
    return NamedFloatingVector(QString(fem_->nodeVectorName(i).c_str()), fem_->nodeVector(i));
}

unsigned long QZScriptEngine::getElementValuesSize() const
{
    return 0;
}

NamedFloatingVector QZScriptEngine::getElementValues(const unsigned long &i)
{
    return NamedFloatingVector();
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

    if (context->argumentCount() == 6 || context->argumentCount() == 7)
    {
        QString typeError = QObject::tr("Quads2D(xCount: Integer, yCount: Integer, origin: Point2D, width: Floating, height: Floating, function: Function[, points: Array]): argument type error (%1).");
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
        return engine->newQObject(new QQuadrilateralMesh2D(xCount, yCount, origin->x(), origin->y(), width, height, func, pointList), QScriptEngine::ScriptOwnership);
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

QScriptValue QZScriptEngine::createQuadrilateralMesh3D(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() == 4)
    {
        QString typeError = QObject::tr("SurfaceQuads(rCount: Integer, lCount: Integer, radius: Floating, length: Floating): argument type error (%1).");
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
        return engine->newQObject(new QQuadrilateralMesh3D(rCount, lCount, radius, length), QScriptEngine::ScriptOwnership);
    }
    return context->throwError(QObject::tr("SurfaceQuads(): arguments count error."));
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
    else if (qscriptvalue_cast<QQuadrilateralMesh3D *>(context->argument(0)) != NULL)
    {
        mesh_ = new QuadrilateralMesh3D(qscriptvalue_cast<QQuadrilateralMesh3D *>(context->argument(0)));
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

QScriptValue QZScriptEngine::planeStress(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() > 4)
    {
        QString typeError = QObject::tr("PlaneStress(mesh: Mesh, h: Floating, E: Floating, nu: Floating, ...): argument type error (%1).");
        if (!context->argument(0).isQObject())
            return context->throwError(typeError.arg("mesh"));
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

        ElasticMatrix D(E, nu, true); //!

        if (fem_ != NULL) delete fem_;

        fem_ = new PlaneStressStrain (dynamic_cast<Mesh2D*>(mesh_), //!
                                      h,
                                      D,
                                      conditions);

        fem_->printNodeValuesExtremums();

        return engine->undefinedValue();;
    }
    return context->throwError(QObject::tr("PlaneStress(mesh: Mesh, h: Floating, E: Floating, nu: Floating, boundaryConditionType: Function, boundaryValue: Function, nodalForce: Function, surfaceForce: Function, volumeForce: Function): arguments count error."));
}

QScriptValue QZScriptEngine::planeStrain(QScriptContext *context, QScriptEngine *engine)
{
    if (context->argumentCount() > 4)
    {
        QString typeError = QObject::tr("PlaneStrain(mesh: Mesh, h: Floating, E: Floating, nu: Floating, ...): argument type error (%1).");
        if (!context->argument(0).isQObject())
            return context->throwError(typeError.arg("mesh"));
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

        ElasticMatrix D(E, nu, false); //!

        if (fem_ != NULL) delete fem_;

        fem_ = new PlaneStressStrain (dynamic_cast<Mesh2D*>(mesh_), //!
                                      h,
                                      D,
                                      conditions);

        fem_->printNodeValuesExtremums();

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
        QString typeError = QObject::tr("MindlinPlate(mesh: Mesh, h: Floating, E: Floating, nu: Floating, ...): argument type error (%1).");
        if (!context->argument(0).isQObject())
            return context->throwError(typeError.arg("mesh"));
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

            ElasticMatrix D(E, nu, true); //!

            if (fem_ != NULL) delete fem_;

            fem_ = new MindlinPlateBending (dynamic_cast<Mesh2D*>(mesh_), //!
                                            h,
                                            D,
                                            conditions);
        }
        else if (context->argument(1).isArray() && context->argument(2).isArray() && context->argument(3).isArray())
        {
            QScriptValue h_array = context->argument(1);
            QScriptValue e_array = context->argument(2);
            QScriptValue nu_array = context->argument(1);
            std::vector<double> h;
            std::vector<ElasticMatrix> elasticMatrix;

            if (h_array.property("length").toInteger() != e_array.property("length").toInteger() || h_array.property("length").toInteger() != nu_array.property("length").toInteger())
            {
                return context->throwError(typeError.arg("h, E, nu: all the input arrays must have same length"));
            }

            for (int i = 0; i < h_array.property("length").toInteger(); i++)
            {
                h.push_back(h_array.property(i).toNumber());
                elasticMatrix.push_back(ElasticMatrix(e_array.property(i).toNumber(), nu_array.property(i).toNumber(), true));
            }

            if (fem_ != NULL) delete fem_;

            fem_ = new MindlinPlateBending (dynamic_cast<Mesh2D*>(mesh_), //!
                                            h,
                                            elasticMatrix,
                                            conditions);
        }
        else
        {
            return context->throwError(typeError.arg("h, E, nu: all the input must be scalars or arrays of same length"));
        }

        fem_->printNodeValuesExtremums();

        return engine->undefinedValue();
    }
    return context->throwError(QObject::tr("PlaneStress(mesh: Mesh, h: Floating, E: Floating, nu: Floating, boundaryConditionType: Function, boundaryValue: Function, nodalForce: Function, surfaceForce: Function, volumeForce: Function): arguments count error."));
}


