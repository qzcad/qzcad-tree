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
/**
 * @brief Qt Script реализация операции приблизительного равенства
 * @param ctx Контекст
 * @param eng Двигатель
 * @return Резульатат операции приблизительного равенства
 */
QScriptValue approx(QScriptContext *context, QScriptEngine *engine);
/**
 * @brief "Конструктор" двумерных точек в скриптах
 * @param context Контекст
 * @param engine Двигатель
 * @return Объект для использования в скриптах
 */
QScriptValue createPoint2D(QScriptContext *context, QScriptEngine *engine);
/**
 * @brief "Конструктор" трехмерных точек в скриптах
 * @param context Контекст
 * @param engine Двигатель
 * @return Объект для использованиия в скриптах
 */
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
