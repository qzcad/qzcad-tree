#include "qtscriptforcecondition3d.h"
#include "iostream"

QtScriptForceCondition3D::QtScriptForceCondition3D(QString condition, QString u, QString v, QString w, ForceType fType)
{
    x_ = y_ = z_ = 0.0;

    setCondition(condition, u, v, w);

    setForceType(fType);
}

bool QtScriptForceCondition3D::setCondition(QString condition, QString u, QString v, QString w)
{
    QString func = "function isApplied(x, y, z){\n";
    func += "return " + condition + ";\n}\n";
    func += "function u(x, y, z){\n";
    func += "return " + u + ";\n}\n";
    func += "function v(x, y, z){\n";
    func += "return " + v + ";\n}\n";
    func += "function w(x, y, z){\n";
    func += "return " + w + ";\n}\n";
    condition_ = condition;
    u_ = u;
    v_ = v;
    w_ = w;
    QScriptSyntaxCheckResult chk = QScriptEngine::checkSyntax(func);
    if (chk.state() != QScriptSyntaxCheckResult::Valid)
    {
        std::cout << "Ошибка анализа граничного условия: " << chk.errorMessage().toStdString() << std::endl;
        return false;
    }
    object_ = engine_.evaluate(func);
    if (object_.isError())return false;
    return true;
}

bool QtScriptForceCondition3D::isU()
{
    return true;
}

double QtScriptForceCondition3D::u()
{
    QScriptValue func = engine_.globalObject().property("u");
    QScriptValue result = func.call(object_, QScriptValueList() << x_ << y_ << z_);
    return result.toNumber();
}

bool QtScriptForceCondition3D::isV()
{
    return true;
}

double QtScriptForceCondition3D::v()
{
    QScriptValue func = engine_.globalObject().property("v");
    QScriptValue result = func.call(object_, QScriptValueList() << x_ << y_ << z_);
    return result.toNumber();
}

bool QtScriptForceCondition3D::isW()
{
    return true;
}

double QtScriptForceCondition3D::w()
{
    QScriptValue func = engine_.globalObject().property("w");
    QScriptValue result = func.call(object_, QScriptValueList() << x_ << y_ << z_);
    return result.toNumber();
}

bool QtScriptForceCondition3D::isApplied(msh::PointPointer point)
{
    QScriptValue func = engine_.globalObject().property("isApplied");
    QScriptValue result = func.call(object_, QScriptValueList() << point->x() << point->y() << point->z());
    x_ = point->x();
    y_ = point->y();
    z_ = point->z();
    return result.toBool();
}
