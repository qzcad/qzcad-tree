/**
  * @author Сергей Чопоров
  * @date 06/08/2015
  * @version 1.0.0
  * @copyright Copyright 2015 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef QTRIANGLEMESH3D_H
#define QTRIANGLEMESH3D_H

#include <QObject>
#include <qmetatype.h>
#include "trianglemesh3d.h"

using namespace msh;
/**
 * @brief Прокси-класс для создания сеток в скриптах моделей
 * @see QObject, TriangleMesh3D
 */
class QTriangleMesh3D : public QObject, public TriangleMesh3D
{
    Q_OBJECT
public:
    /**
     * @brief Конструктор по умолчанию
     * @param parent Указатель на родительский объект
     */
    explicit QTriangleMesh3D(QObject *parent = 0);
    /**
     * @brief Конструктор копирования
     * @param mesh Экземпляр объекта для копирования
     */
    explicit QTriangleMesh3D(const QTriangleMesh3D &mesh);
    /**
     * @brief Конструктор для создания равномерной сетки в цилиндрических координатах
     * @param rCount Количество элементов вдоль радиуса
     * @param lCount Количество элементов вдоль образующей
     * @param radius Радиус
     * @param length Длина образующей
     * @param parent Указатель на родительский объект
     */
    explicit QTriangleMesh3D(const UInteger &rCount, const UInteger &lCount, const double &radius, const double &length, QObject *parent = 0);
    /**
     * @brief Конструктор для создания равномерной сетки в конических координатах
     * @param rCount Количество элементов вдоль радиуса
     * @param lCount Количество элементов вдоль образующей
     * @param bottom_radius Нижний радиус
     * @param top_radius Верхний радиус
     * @param length Длина образующей
     * @param parent Указатель на родительский объект
     */
    explicit QTriangleMesh3D(const UInteger &rCount, const UInteger &lCount, const double &bottom_radius, const double &top_radius, const double &length, QObject *parent = 0);
    /**
     * @brief Метод формирует строку с информацией о сетке
     * @return Строка с информацией о сетке
     * Метод доступен в скриптах
     */
    Q_INVOKABLE QString toString() const;
signals:

public slots:

};
// Регистрация класса для использования в скриптах моделей
Q_DECLARE_METATYPE(QTriangleMesh3D)
Q_DECLARE_METATYPE(QTriangleMesh3D *)
#endif // QTRIANGLEMESH3D_H
