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

#include "fem.h"

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
    static double epsilon();
    /**
     * @brief Установить значение точности (машинного нуля) вычислений (минимальное растояние между двумя неравными числами)
     * @param epsilon Новое значение точности (машинного нуля) вычислений (минимальное растояние между двумя неравными числами)
     */
    static void setEpsilon(double epsilon);
    /**
     * @brief Получить указатель на созданую сетку
     * @return Указатель на сетку
     */
    static msh::Mesh *mesh();
    /**
     * @brief Установить указатель текущей сетки
     * @param mesh Указатель на сетку
     */
    static void setMesh(msh::Mesh *mesh);
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
     * @brief "Конструктор" точек в скриптах
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Объект для использования в скриптах
     */
    static QScriptValue createPoint(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief "Конструктор" сеток четырехугольных элементов на плоскости
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Объект для использования в скриптах
     */
    static QScriptValue createQuadrilateralMesh2D(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief "Конструктор" сеток двумерных балок
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Объект для использования в скриптах
     */
    static QScriptValue createSegmentMesh2D(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief "Конструктор" сеток треугольных элементов на плоскости
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Объект для использования в скриптах
     */
    static QScriptValue createTriangleMesh2D(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief "Конструктор" сеток треугольных при помощи триангуляции Делоне
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Объект для использования в скриптах
     */
    static QScriptValue createDelaunay(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief "Конструктор" сеток треугольных при помощи триангуляции Делоне с использование сглаживани Рапперта
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Объект для использования в скриптах
     */
    static QScriptValue createRuppert(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief "Конструктор" сеток четырехугольных элементов в пространстве (оболочка)
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Объект для использования в скриптах
     */
    static QScriptValue createQuadrilateralMesh3D(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief "Конструктор" сеток четырехугольных элементов в цилиндрических координатах
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Объект для использования в скриптах
     */
    static QScriptValue createCylinderQuads(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief "Конструктор" сеток четырехугольных элементов в конических координатах (усеченный конус)
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Объект для использования в скриптах
     */
    static QScriptValue createConeQuads(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief "Конструктор" сеток треугольных элементов в пространстве (оболочка)
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Объект для использования в скриптах
     */
    static QScriptValue createTriangleMesh3D(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief "Конструктор" сеток треугольных элементов в цилиндрических координатах
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Объект для использования в скриптах
     */
    static QScriptValue createCylinderTriangles(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief "Конструктор" сеток треугольных элементов в конических координатах
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Объект для использования в скриптах
     */
    static QScriptValue createConeTriangles(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief Qt Script реализация операции печати на консоль через std::cout
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return engine->undefinedValue()
     */
    static QScriptValue printStd(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief Установить результирующую сетку
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return engine->undefinedValue()
     */
    static QScriptValue setMesh(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief Получить копию текущей сетки в перменной
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Соответсвующий указатель
     */
    static QScriptValue currentMesh(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief Функция суммирования с произвольным числом аргументов
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Если аргументы корректны, то их сумму, иначе engine->undefinedValue()
     */
    static QScriptValue sum(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief Конъюнкция произвольного числа аргументов
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Если аргументы корректны (тип Number), то значение конъюнкции для них, иначе engine->undefinedValue()
     */
    static QScriptValue con(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief Дизъюнкция произвольного числа аргументов
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Если аргументы корректны (тип Number), то значение дизъюнкции для них, иначе engine->undefinedValue()
     */
    static QScriptValue dis(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief Логическая разность произвольного числа аргументов: первый минус остальные
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Если аргументы корректны (тип Number), то значение дизъюнкции для них, иначе engine->undefinedValue()
     */
    static QScriptValue diff(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief Область, ограниченная окружностью
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Если аргументы корректны (тип Number), то значение соответствующей функциий, иначе engine->undefinedValue()
     */
    static QScriptValue circle(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief Область, ограниченная эллипсом
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Если аргументы корректны (тип Number), то значение соответствующей функциий, иначе engine->undefinedValue()
     */
    static QScriptValue ellipse(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief Полоса, перепендикулярная оси, симметричная относительно начала отсчета
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Если аргументы корректны (тип Number), то значение соответствующей функциий, иначе engine->undefinedValue()
     */
    static QScriptValue band(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief Область находящаяся "выше" прямой, определенной точками (x1, y1), (x2, y2)
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Если аргументы корректны (тип Number), то значение соответствующей функциий, иначе engine->undefinedValue()
     */
    static QScriptValue line(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief Функция области, ограниченной прямоугольником с центром в начале координат
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Если аргументы корректны (тип Number), то значение соответствующей функциий, иначе engine->undefinedValue()
     */
    static QScriptValue rectangle(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief МКЭ: плоское напряжение
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Если аргументы корректны, то будет найдено решение задачи, иначе engine->undefinedValue()
     */
    static QScriptValue planeStress(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief МКЭ: плоская деформация
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Если аргументы корректны, то будет найдено решение задачи, иначе engine->undefinedValue()
     */
    static QScriptValue planeStrain(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief Функция регистрации граничного условия
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Если аргументы корректны, то зарегистрированное граничное условие, иначе engine->undefinedValue()
     */
    static QScriptValue createBoundaryCondition(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief Функция регистрации узловой нагрузки
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Если аргументы корректны, то зарегистрированное граничное условие, иначе engine->undefinedValue()
     */
    static QScriptValue createNodalForce(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief Функция регистрации поверхностной нагрузки
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Если аргументы корректны, то зарегистрированное граничное условие, иначе engine->undefinedValue()
     */
    static QScriptValue createSurfaceForce(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief Функция регистрации объемной нагрузки
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Если аргументы корректны, то зарегистрированное граничное условие, иначе engine->undefinedValue()
     */
    static QScriptValue createVolumeForce(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief МКЭ: Прогиб пластинки (теория Миндлина)
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Если аргументы корректны, то будет найдено решение задачи, иначе engine->undefinedValue()
     */
    static QScriptValue mindlinPlate(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief МКЭ: Деформация оболочки (теория Миндлина)
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Если аргументы корректны, то будет найдено решение задачи, иначе engine->undefinedValue()
     */
    static QScriptValue mindlinShell(QScriptContext *context, QScriptEngine *engine);
    /**
     * @brief Отчет о значениях
     * @param context Контекст скрипта
     * @param engine Двигатель скрипта
     * @return Если аргументы корректны, то будет выведен в стандартный вывод отчет
     */
    static QScriptValue reportValues(QScriptContext *context, QScriptEngine *engine);
private:
    static double epsilon_; //!< Точность численных операций
    static msh::Mesh *mesh_; //!< Сетка, построенная в результате интерпретации скрипта
    static Fem *fem_; //!< Результаты МКЭ
};

#endif // QZSCRIPTENGINE_H
