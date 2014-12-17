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

#include "mesh.h"

#include "namedfloatingvector.h"

/**
 * @brief Интерпретатор скриптов модели
 * @see QScriptEngine
 */
class QZScriptEngine : public QScriptEngine
{
    Q_OBJECT
public:
    /**
     * @brief Конструктор по умолчанию
     * @param parent Указатель на родительский объект
     */
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
    /**
     * @brief Получить указатель на созданую сетку
     * @return Указатель на сетку
     */
    msh::Mesh *mesh();
    /**
     * @brief Получить количество векторов со значениями, определенными в узлах
     * @return Количество векторов со значениями, определенными в узлах
     */
    unsigned long getNodeValuesSize() const;
    /**
     * @brief Получить вектор, значений опеределенных в узлах
     * @param i Номер вектора
     * @return Вектор, значений определенных в узлах
     */
    NamedFloatingVector &getNodeValues(const unsigned long &i);
    /**
     * @brief Получить количество векторов со значениями, определенными на элементах
     * @return Количество векторов со значениями, определенными на элементах
     */
    unsigned long getElementValuesSize() const;
    /**
     * @brief Получить вектор, значений опеределенных на элементах
     * @param i Номер вектора
     * @return Вектор, значений определенных на элементах
     */
    NamedFloatingVector &getElementValues(const unsigned long &i);
signals:

public slots:
private:
    /**
     * @brief Функция для печати на екран сообщения об версии интерпретатора
     * @param context Контекст скрипта
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
    msh::Mesh *mesh_; //!< Сетка, построенная в результате интерпретации скрипта
    std::vector<NamedFloatingVector> nodeValues_; //!< Массив векторов значений, определенных в узле сетки
    std::vector<NamedFloatingVector> elementValues_; //!< Массив векторов значений, определенных в узле сетки
};

#endif // QZSCRIPTENGINE_H
