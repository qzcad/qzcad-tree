/**
  * @author Сергей Чопоров
  * @date 12/12/2014
  * @version 1.0.0
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef QPOINT3D_H
#define QPOINT3D_H

#include <QObject>
#include <qmetatype.h>
#include "point3d.h"

using namespace msh;
/**
 * @brief Прокси-класс для декларирования и использования точек трехмерного пространства в скриптах моделей
 * @see QObject, Point3D
 */
class QPoint3D : public QObject, public Point3D
{
    Q_OBJECT
public:
    /**
     * @brief Конструктор по умолчанию
     * @param parent Указатель на родительский объект
     */
    explicit QPoint3D(QObject *parent = 0);
    /**
     * @brief Конструктор для инициализации координат
     * @param x Абсцисса
     * @param y Ордината
     * @param z Аппликата
     * @param parent Указатель на родительский объект
     */
    QPoint3D(const double &x, const double &y, const double &z, QObject *parent = 0);
    /**
     * @brief Конструктор копирования
     * @param point Экземпляр объекта для копирования
     */
    QPoint3D(const QPoint3D &point);
    /**
     * @brief Преобразование в строку вида "(x; y; z)"
     * @return Строка с координатами точки
     */
    Q_INVOKABLE QString toString() const;
    /**
     * @brief Векторное произведение векторов, определенных текущей точкой и заданной
     * @param point Координаты ветора, с которым необходимо найти векторное произведение
     * @return Значение векторного произведения
     */
    Q_INVOKABLE Point3D product(const QPoint3D &point) const;
    /// Объявление полей, доступных в скриптах
    Q_PROPERTY(double x READ x WRITE setX) //!< Абсцисса
    Q_PROPERTY(double y READ y WRITE setY) //!< Ордината
    Q_PROPERTY(double z READ z WRITE setZ) //!< Аппликата
signals:

public slots:

};
// Регистрация класса для использования в скриптах моделей
Q_DECLARE_METATYPE(QPoint3D)
Q_DECLARE_METATYPE(QPoint3D*)

#endif // QPOINT3D_H
