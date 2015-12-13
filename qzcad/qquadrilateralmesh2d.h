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
     * @brief Конструктор копирования
     * @param mesh Экземпляр объекта для копирования
     */
    QQuadrilateralMesh2D(const QQuadrilateralMesh2D &mesh);
    /**
     * @brief Конструктор создает копию объекта, переданного по указателю
     * @param mesh Указатель на объект для копирования
     * @param parent Указатель на родительский объект
     */
    QQuadrilateralMesh2D(QuadrilateralMesh2D *mesh, QObject *parent = 0);
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
