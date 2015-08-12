#include "qfemcondition.h"

QFemCondition::QFemCondition(QObject *parent) :
    QObject(parent)
{
//    value_ = 0.0;
    type_ = INITIAL_VALUE;
    direction_ = ALL;
//    isConst_ = true;
}

//QFemCondition::QFemCondition(int type, int direction, double value, const QScriptValue &conditionFunction, QObject *parent) :
//    QObject(parent)
//{
//    type_ = type;
//    direction_ = direction;
//    conditionFunction_ = conditionFunction;
//    value_ = value;
//    isConst_ = true;
//}

QFemCondition::QFemCondition(FemConditionType type, FemDirection direction, const QScriptValue &conditionFunction, const QScriptValue &valueFunction, QObject *parent) :
    QObject(parent)
{
    type_ = type;
    direction_ = direction;
    conditionFunction_ = conditionFunction;
    valueFunction_ = valueFunction;
//    isConst_ = false;
}

QFemCondition::QFemCondition(const QFemCondition &qfc) :
    QObject(qfc.parent())
{
    type_ = qfc.type_;
    direction_ = qfc.direction_;
//    value_ = qfc.value_;
//    isConst_ = qfc.isConst_;
    conditionFunction_ = qfc.conditionFunction_;
    valueFunction_ = qfc.valueFunction_;
}

FemCondition::FemConditionType QFemCondition::type() const
{
    return type_;
}

FemCondition::FemDirection QFemCondition::direction() const
{
    return direction_;
}

bool QFemCondition::isApplied(msh::PointPointer point)
{
    if (conditionFunction_.isFunction())
    {
        QScriptValueList args;
        args << point->x() << point->y() << point->z();
        return conditionFunction_.call(QScriptValue(), args).toBool();
    }
    return conditionFunction_.toBool();
}

double QFemCondition::value(msh::PointPointer point)
{
//    if (isConst_)
//        return value_;

//    QScriptValueList args;
//    args << point->x() << point->y() << point->z();
//    return valueFunction_.call(QScriptValue(), args).toNumber();
    if (valueFunction_.isFunction())
    {
        QScriptValueList args;
        args << point->x() << point->y() << point->z();
        return valueFunction_.call(QScriptValue(), args).toNumber();
    }
    return valueFunction_.toNumber();
}

QString QFemCondition::toString() const
{
    return QString("type = %1;\tdirection = %2").arg(type_).arg(direction_);
}
