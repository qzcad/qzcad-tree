/**
  * @author Сергей Чопоров
  * @date 13/04/2015
  * @version 1.0.0
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef QTRIANGLEMESH2D_H
#define QTRIANGLEMESH2D_H

#include <QObject>
#include <qmetatype.h>
#include "trianglemesh2d.h"

using namespace msh;

/**
 * @brief Прокси-класс для создания сеток в скриптах моделей
 * @see QObject, QuadrilateralMesh2D
 */
class QTriangleMesh2D : public QObject, public TriangleMesh2D
{
    Q_OBJECT
public:
    /**
     * @brief Конструктор по умолчанию
     * @param parent Указатель на родительский объект
     */
    explicit QTriangleMesh2D(QObject *parent = 0);
    /**
     * @brief Конструктор создает равномерную структурированную секту в прямоугольной области
     * @param xCount Количество узлов вдоль оси абсцисс
     * @param yCount Количество узлов вдоль оси ординат
     * @param xMin Абсцисса нижнего левого угла прямоугольной области
     * @param yMin Ордината нижнего левого угла прямоугольной области
     * @param width Ширина прямоугольной области
     * @param height Высота прямоугольной области
     * @param parent Указатель на родительский объект
     */
    explicit QTriangleMesh2D(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, QObject *parent = 0);
    /**
     * @brief Конструктор копирования
     * @param qmesh Экземпляр объекта для копирования
     */
    explicit QTriangleMesh2D(const QTriangleMesh2D &qmesh);
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
Q_DECLARE_METATYPE(QTriangleMesh2D)
Q_DECLARE_METATYPE(QTriangleMesh2D*)
#endif // QTRIANGLEMESH2D_H
