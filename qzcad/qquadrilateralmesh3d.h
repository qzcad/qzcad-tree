/**
  * @author Сергей Чопоров
  * @date 06/08/2015
  * @version 1.0.0
  * @copyright Copyright 2015 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef QQUADRILATERALMESH3D_H
#define QQUADRILATERALMESH3D_H

#include <QObject>
#include <qmetatype.h>
#include "quadrilateralmesh3d.h"

using namespace msh;
/**
 * @brief Прокси-класс для создания сеток в скриптах моделей
 * @see QObject, QuadrilateralMesh3D
 */
class QQuadrilateralMesh3D : public QObject, public QuadrilateralMesh3D
{
    Q_OBJECT
public:
    /**
     * @brief Конструктор по умолчанию
     * @param parent Указатель на родительский объект
     */
    explicit QQuadrilateralMesh3D(QObject *parent = 0);
    /**
     * @brief Конструктор копирования
     * @param mesh Экземпляр объекта для копирования
     */
    QQuadrilateralMesh3D(const QQuadrilateralMesh3D &mesh);
    /**
     * @brief Конструктор для создания равномерной сетки в цилиндрических координатах
     * @param rCount Количество элементов вдоль радиуса
     * @param lCount Количество элементов вдоль образующей
     * @param radius Радиус
     * @param length Длина образующей
     * @param parent Указатель на родительский объект
     */
    QQuadrilateralMesh3D(const UInteger &rCount, const UInteger &lCount, const double &radius, const double &length, QObject *parent = 0);
    /**
     * @brief Конструктор для создания равномерной сетки в конических координатах
     * @param rCount Количество элементов вдоль радиуса
     * @param lCount Количество элементов вдоль образующей
     * @param bottom_radius Нижний радиус
     * @param top_radius Верхний радиус
     * @param length Длина образующей
     * @param parent Указатель на родительский объект
     */
    QQuadrilateralMesh3D(const UInteger &rCount, const UInteger &lCount, const double &bottom_radius, const double &top_radius, const double &length, QObject *parent = 0);
    /**
     * @brief Конструктор создает копию объекта, переданного по указателю
     * @param mesh Указатель на объект для копирования
     * @param parent Указатель на родительский объект
     */
    QQuadrilateralMesh3D(QuadrilateralMesh3D *mesh, QObject *parent = 0);
    /**
     * @brief Метод формирует строку с информацией о сетке
     * @return Строка с информацией о сетке
     * Метод доступен в скриптах
     */
    Q_INVOKABLE QString toString() const;
    /**
     * @brief Операция добавления сетки к существующей на основе элементарного объединения
     * @param mesh Указатель на сетку для добавления
     */
    Q_INVOKABLE void add(const QQuadrilateralMesh3D *mesh);
    /**
     * @brief Операция перемещения сетки на заданный радиус вектор
     * @param x Абсцисса радиус вектора перемещения
     * @param y Ордината радиус вектора перемещения
     * @param y Аппликата радиус вектора перемещения
     */
    Q_INVOKABLE void translate(const double &x, const double &y, const double &z);
signals:

public slots:

};
// Регистрация класса для использования в скриптах моделей
Q_DECLARE_METATYPE(QQuadrilateralMesh3D)
Q_DECLARE_METATYPE(QQuadrilateralMesh3D *)
#endif // QQUADRILATERALMESH3D_H
