/**
  * @author Сергей Чопоров
  * @date 07/12/2014
  * @version 1.0.0
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef QZSCRIPTENGINE_H
#define QZSCRIPTENGINE_H

#include <QScriptEngine>
/**
 * @brief Интерпретатор скриптов модели
 * @see QScriptEngine
 */
class QZScriptEngine : public QScriptEngine
{
    Q_OBJECT
public:
    explicit QZScriptEngine(QObject *parent = 0);
    /**
     * @brief Получить значение точности (машинного нуля) вычислений (минимальное растояние между двумя неравными числами)
     * @return Текущее значение точности
     */
    double epsilon() const;
    /**
     * @brief Установить значение точности (машинного нуля) вычислений (минимальное растояние между двумя неравными числами)
     * @param epsilon Новое значение точности (машинного нуля) вычислений (минимальное растояние между двумя неравными числами)
     */
    void setEpsilon(double epsilon);

signals:

public slots:
private:
    /**
     * @brief Функция для печати на екран сообщения об версии интерпретатора
     * @param context Контекст скрипт
     * @param engine Двигатель скрипта
     * @return engine->undefinedValue()
     */
    static QScriptValue about(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief Qt Script реализация операции приблизительного равенства
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Резульатат операции приблизительного равенства
     */
    static QScriptValue approx(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief "Конструктор" двумерных точек в скриптах
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Объект для использования в скриптах
     */
    static QScriptValue createPoint2D(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief "Конструктор" трехмерных точек в скриптах
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Объект для использованиия в скриптах
     */
    static QScriptValue createPoint3D(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief Qt Script реализация операции печати на консоль через std::cout
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return engine->undefinedValue()
     */
    static QScriptValue printStd(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief Функция суммирования с произвольным числом аргументов
     * @param context Конетекст скрипта
     * @param engine Двигатель скрипта
     * @return Если аргументы корректны, то их сумму, иначе engine->undefinedValue()
     */
    static QScriptValue sum(QScriptContext *context, QScriptEngine *engine);
private:
    static double epsilon_; //!< Точность численных операций
};

#endif // QZSCRIPTENGINE_H
