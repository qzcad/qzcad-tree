#ifndef QFEMCONDITION_H
#define QFEMCONDITION_H

#include <QObject>
#include <QScriptEngine>
#include <qmetatype.h>
#include "femcondition.h"
#include "pointpointer.h"

class QFemCondition : public QObject, public FemCondition
{
    Q_OBJECT
public:
    explicit QFemCondition(QObject *parent = 0);
//    QFemCondition(int type, int direction, double value, const QScriptValue &conditionFunction, QObject *parent = 0);
    QFemCondition(FemConditionType type, int direction, const QScriptValue &conditionFunction, const QScriptValue &valueFunction, QObject *parent = 0);
    QFemCondition(const QFemCondition &qfc);
    FemConditionType type() const;
    int direction() const;
    bool isApplied(msh::PointPointer point);
    double value(msh::PointPointer point);
    Q_INVOKABLE QString toString() const;
signals:

public slots:
private:
    QScriptValue conditionFunction_; //!< Функция-условие для отбора точек
    QScriptValue valueFunction_; //!< Функция - значение в точке
//    double value_; //!< Начальное значение для случая, когда функция принимает константное значение по всей области определения
    FemConditionType type_;
    int direction_;
//    bool isConst_;
};
// Регистрация класса для использования в скриптах моделей
Q_DECLARE_METATYPE(QFemCondition)
Q_DECLARE_METATYPE(QFemCondition*)
#endif // QFEMCONDITION_H
