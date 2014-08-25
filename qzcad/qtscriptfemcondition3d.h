/**
  * @author Сергей Чопоров
  * @date 24/08/2014
  * @version 1.0.0
  **/
#ifndef QTSCRIPTFEMCONDITION3D_H
#define QTSCRIPTFEMCONDITION3D_H
#include <QtScript/QScriptEngine>
#include "femcondition3d.h"
/**
 * @brief Граничные условия на базе QtScript
 * @see QtScript
 */
class QtScriptFemCondition3D : public FEMCondition3D
{
public:
    QtScriptFemCondition3D(QString condition, bool isU, double u, bool isV, double v, bool isW, double w);
    QtScriptFemCondition3D(const QtScriptFemCondition3D &femCondition);
    bool setCondition(QString script);
    virtual bool isU();
    void setIsU(bool isU);

    virtual double u();
    void setU(double u);

    virtual bool isV();
    void setIsV(bool isV);

    virtual double v();
    void setV(double v);

    virtual bool isW();
    void setIsW(bool isW);

    virtual double w();
    void setW(double w);
    virtual bool isApplied(msh::PointPointer point);
    virtual ~QtScriptFemCondition3D(){}
private:
    QScriptEngine engine_;
    QScriptValue object_;
    bool isU_;
    double u_;
    bool isV_;
    double v_;
    bool isW_;
    double w_;
    QString script_;
};

#endif // QTSCRIPTFEMCONDITION3D_H
