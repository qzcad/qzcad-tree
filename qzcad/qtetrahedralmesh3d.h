/**
  * @author Сергей Чопоров
  * @date 22/03/2017
  * @version 1.0.0
  * @copyright Copyright 2017 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef QTETRAHEDRALMESH3D_H
#define QTETRAHEDRALMESH3D_H

#include <QObject>
#include <qmetatype.h>
#include "tetrahedralmesh3d.h"

using namespace msh;

/**
 * @brief Прокси-класс для создания сеток в скриптах моделей
 * @see QObject, TetrahedralMesh3D
 */
class QTetrahedralMesh3D : public QObject, public TetrahedralMesh3D
{
    Q_OBJECT
public:
    /**
     * @brief Конструктор по умолчанию
     * @param parent Указатель на родительский объект
     */
    explicit QTetrahedralMesh3D(QObject *parent = 0);
    /**
     * @brief Конструктор копирования
     * @param qmesh Экземпляр объекта для копирования
     */
    QTetrahedralMesh3D(const QTetrahedralMesh3D &qmesh);
    /**
     * @brief Конструктор создает копию объекта, переданного по указателю
     * @param mesh Указатель на объект для копирования
     * @param parent Указатель на родительский объект
     */
    QTetrahedralMesh3D(TetrahedralMesh3D *mesh, QObject *parent = 0);
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
Q_DECLARE_METATYPE(QTetrahedralMesh3D)
Q_DECLARE_METATYPE(QTetrahedralMesh3D*)
#endif // QTETRAHEDRALMESH3D_H
