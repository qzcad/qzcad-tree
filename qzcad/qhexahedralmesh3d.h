/**
  * @author Сергей Чопоров
  * @date 02/10/2017
  * @version 1.0.0
  * @copyright Copyright 2017 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef QHEXAHEDRALMESH3D_H
#define QHEXAHEDRALMESH3D_H

#include <QObject>
#include <qmetatype.h>
#include "hexahedralmesh3d.h"

using namespace msh;

class QHexahedralMesh3D : public QObject,  public HexahedralMesh3D
{
    Q_OBJECT
public:
    /**
     * @brief Конструктор по умолчанию
     * @param parent Указатель на родительский объект
     */
    explicit QHexahedralMesh3D(QObject *parent = 0);
    /**
     * @brief Конструктор копирования
     * @param qmesh Экземпляр объекта для копирования
     */
    QHexahedralMesh3D(const QHexahedralMesh3D &qmesh);
    /**
     * @brief Конструктор создает копию объекта, переданного по указателю
     * @param mesh Указатель на объект для копирования
     * @param parent Указатель на родительский объект
     */
    QHexahedralMesh3D(HexahedralMesh3D *mesh, QObject *parent = 0);
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
Q_DECLARE_METATYPE(QHexahedralMesh3D)
Q_DECLARE_METATYPE(QHexahedralMesh3D*)
#endif // QHEXAHEDRALMESH3D_H
