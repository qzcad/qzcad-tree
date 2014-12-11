/**
  * @author Сергей Чопоров
  * @date 11/12/2014
  * @version 1.0.0
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef QPOINT2D_H
#define QPOINT2D_H

#include <QObject>
#include <qmetatype.h>
#include "point2d.h"

using namespace msh;
/**
 * @brief Прокси-класс для декларирования и использования точек на плоскости в скриптах моделей
 */
class QPoint2D : public QObject, public Point2D
{
    Q_OBJECT
public:
    /**
     * @brief Конструктор по умолчанию
     * @param parent Указатель на родительский объект
     */
    explicit QPoint2D(QObject *parent = 0);
    /**
     * @brief Конструктор для инициализации координат
     * @param x Абсцисса
     * @param y Ордината
     * @param parent Указатель на родительский объект
     */
    explicit QPoint2D(double x, double y, QObject *parent = 0);
    /**
     * @brief Конструктор копирования
     * @param point Экземпляр для копирования
     */
    explicit QPoint2D(const QPoint2D &point);
    /**
     * @brief Преобразовать в строку вида "(x; y)"
     * @return Строка с координатами точки
     */
    Q_INVOKABLE QString toString() const;
    /// Объявление полей, доступных в скриптах
    Q_PROPERTY(double x READ x WRITE setX) //!< Абсцисса
    Q_PROPERTY(double y READ y WRITE setY) //!< Ордината
signals:

public slots:
};
// Регистрация класса для использования в скриптах моделей
Q_DECLARE_METATYPE(QPoint2D)
Q_DECLARE_METATYPE(QPoint2D*)

#endif // QPOINT2D_H
