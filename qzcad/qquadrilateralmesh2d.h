/**
  * @author Сергей Чопоров
  * @date 09/12/2014
  * @version 1.0.0
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef QQUADRILATERALMESH2D_H
#define QQUADRILATERALMESH2D_H

#include <QObject>
#include <qmetatype.h>
#include "quadrilateralmesh2d.h"

using namespace msh;
/**
 * @brief Прокси-класс для создания сеток в скриптах моделей
 * @see QObject, QuadrilateralMesh2D
 */
class QQuadrilateralMesh2D : public QObject, public QuadrilateralMesh2D
{
    Q_OBJECT
public:
    /**
     * @brief Конструктор по умолчанию
     * @param parent Указатель на родительский объект
     */
    explicit QQuadrilateralMesh2D(QObject *parent = 0);
    /**
     * @brief Конструктор создает равномерную структурированную секту в прямоугольной областе
     * @param xCount Количество узлов вдоль оси абсцисс
     * @param yCount Количество узлов вдоль оси ординат
     * @param xMin Абсцисса нижнего левого угла прямоугольной области
     * @param yMin Ордината нижнего левого угла прямоугольной области
     * @param width Ширина прямоугольной области
     * @param height Высота прямоугольной области
     * @param parent Указатель на родительский объект
     */
    explicit QQuadrilateralMesh2D(const UInteger &xCount, const UInteger &yCount, const double &xMin, const double &yMin, const double &width, const double &height, QObject *parent = 0);
    /**
     * @brief Конструктор создает равномерную структурированную сетку для выпуклой четырехугольной области
     * @param xCount Количество узлов по первому направлению (ребра 0-1 и 2-3)
     * @param yCount Количество узлов по второму направлению (ребра 1-2 и 3-0)
     * @param v0 Координаты узла 0
     * @param v1 Координаты узла 1
     * @param v2 Координаты узла 2
     * @param v3 Координаты узла 3
     * @param parent Указатель на родительский объект
     */
    explicit QQuadrilateralMesh2D(const UInteger &xCount, const UInteger &yCount, const Point2D &v0, const Point2D &v1, const Point2D &v2, const Point2D &v3, QObject *parent = 0);
    /**
     * @brief Конструктор создает сетку для треугольной области
     * @param count Количество узлов на сторону треугольника (должно быть четным)
     * @param v0 Координаты узла 0
     * @param v1 Координаты узла 1
     * @param v2 Координаты узла 2
     * @param parent Указатель на родительский объект
     */
    explicit QQuadrilateralMesh2D(const UInteger &count, const Point2D &v0, const Point2D &v1, const Point2D &v2, QObject *parent = 0);
    /**
     * @brief Конструктор копирования
     * @param mesh Экземпляр объекта для копирования
     */
    explicit QQuadrilateralMesh2D(const QQuadrilateralMesh2D &mesh);
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
Q_DECLARE_METATYPE(QQuadrilateralMesh2D)
Q_DECLARE_METATYPE(QQuadrilateralMesh2D *)
#endif // QQUADRILATERALMESH2D_H
