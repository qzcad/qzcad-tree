#include "qtscriptfunctions.h"
#include "qzscriptengine.h"

QZScriptEngine::QZScriptEngine(QObject *parent) :
    QScriptEngine(parent)
{
    // Операция приблизительного равентсва
    QScriptValue qsApprox = newFunction(approx);
    globalObject().setProperty("approx", qsApprox);
    // Переопределение печати
    QScriptValue qsPrint = newFunction(printStd);
    globalObject().setProperty("print", qsPrint);
    // Функция суммирования
    QScriptValue qsSum = newFunction(sum);
    globalObject().setProperty("sum", qsSum);
    // Точка на плоскости
    QScriptValue qsCreatePoint2D = newFunction(createPoint2D);
    globalObject().setProperty("Point2D", qsCreatePoint2D);
    qScriptRegisterMetaType(this, toScriptValuePoint2D, fromScriptValuePoint2D);
    // Точка в пространстве
    QScriptValue qsCreatePoint3D = newFunction(createPoint3D);
    globalObject().setProperty("Point3D", qsCreatePoint3D);
    qScriptRegisterMetaType(this, toScriptValuePoint3D, fromScriptValuePoint3D);
}
