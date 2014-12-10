/**
  * @author Сергей Чопоров
  * @date 07/12/2014
  * @version 1.0.1
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef QTSCRIPTFUNCTIONS_H
#define QTSCRIPTFUNCTIONS_H
#include <QtScript>
#include "point2d.h"
#include "point3d.h"
/**
 * @brief Qt Script реализация операции приблизительного равенства
 * @param ctx Контекст
 * @param eng Двигатель
 * @return Резульатат операции приблизительного равенства
 */
QScriptValue approx(QScriptContext *context, QScriptEngine *engine);
// Point2D
Q_DECLARE_METATYPE(msh::Point2D)
Q_DECLARE_METATYPE(msh::Point2D*)
QScriptValue toScriptValuePoint2D(QScriptEngine *engine, const msh::Point2D &point);
void fromScriptValuePoint2D(const QScriptValue &value, msh::Point2D &point);
QScriptValue createPoint2D(QScriptContext *context, QScriptEngine *engine);
// Point3D
Q_DECLARE_METATYPE(msh::Point3D)
Q_DECLARE_METATYPE(msh::Point3D*)
QScriptValue toScriptValuePoint3D(QScriptEngine *engine, const msh::Point3D &point);
void fromScriptValuePoint3D(const QScriptValue &value, msh::Point3D &point);
QScriptValue createPoint3D(QScriptContext *context, QScriptEngine *engine);
/**
 * @brief Qt Script реализация операции печати на консоль через std::cout
 * @param context Контекст
 * @param engine Двигатель
 * @return eng->undefinedValue()
 */
QScriptValue printStd(QScriptContext *context, QScriptEngine *engine);
QScriptValue sum(QScriptContext *context, QScriptEngine *engine);
#endif // QTSCRIPTFUNCTIONS_H
