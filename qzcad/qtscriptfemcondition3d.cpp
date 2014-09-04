#include "qtscriptfemcondition3d.h"
#include "iostream"
#include "qtscriptfunctions.h"

QtScriptFemCondition3D::QtScriptFemCondition3D(QString condition, bool isU, double u, bool isV, double v, bool isW, double w)
{
    isU_ = isU;
    u_ = u;
    isV_ = isV;
    v_ = v;
    isW_= isW;
    w_ = w;
    script_ = condition;
    setCondition(script_);

    QScriptValue qsApprox = engine_.newFunction(approx);
    engine_.globalObject().setProperty("approx", qsApprox);
}

QtScriptFemCondition3D::QtScriptFemCondition3D(const QtScriptFemCondition3D &femCondition)
{
    isU_ = femCondition.isU_;
    u_ = femCondition.u_;
    isV_ = femCondition.isV_;
    v_ = femCondition.v_;
    isW_ = femCondition.isW_;
    w_ = femCondition.w_;
    script_ = femCondition.script_;
    setCondition(script_);

    QScriptValue qsApprox = engine_.newFunction(approx);
    engine_.globalObject().setProperty("approx", qsApprox);
}

bool QtScriptFemCondition3D::setCondition(QString script)
{
    QString func = "function isApplied(x, y, z){\n";
    func += "return " + script + ";\n}";
    script_ = script;
    QScriptSyntaxCheckResult chk = QScriptEngine::checkSyntax(func);
    if (chk.state() != QScriptSyntaxCheckResult::Valid)
    {
        std::cout << "Ошибка анализа граничного условия: " << chk.errorMessage().toAscii().data() << std::endl;
        return false;
    }
    object_ = engine_.evaluate(func);
    if (object_.isError())return false;
    return true;
}

bool QtScriptFemCondition3D::isU()
{
    return isU_;
}

void QtScriptFemCondition3D::setIsU(bool isU)
{
    isU_ = isU;
}
double QtScriptFemCondition3D::u()
{
    return u_;
}

void QtScriptFemCondition3D::setU(double u)
{
    u_ = u;
}
bool QtScriptFemCondition3D::isV()
{
    return isV_;
}

void QtScriptFemCondition3D::setIsV(bool isV)
{
    isV_ = isV;
}
double QtScriptFemCondition3D::v()
{
    return v_;
}

void QtScriptFemCondition3D::setV(double v)
{
    v_ = v;
}
bool QtScriptFemCondition3D::isW()
{
    return isW_;
}

void QtScriptFemCondition3D::setIsW(bool isW)
{
    isW_ = isW;
}
double QtScriptFemCondition3D::w()
{
    return w_;
}

void QtScriptFemCondition3D::setW(double w)
{
    w_ = w;
}

bool QtScriptFemCondition3D::isApplied(msh::PointPointer point)
{
    QScriptValue func = engine_.globalObject().property("isApplied");
    QScriptValue result = func.call(object_, QScriptValueList() << point->x() << point->y() << point->z());
    return result.toBool();
}






