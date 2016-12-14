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
    QTriangleMesh3D(const QTriangleMesh3D &mesh);
    /**
     * @brief Конструктор создает копию объекта, переданного по указателю
     * @param mesh Указатель на объект для копирования
     * @param parent Указатель на родительский объект
     */
    QTriangleMesh3D(TriangleMesh3D *mesh, QObject *parent = 0);
    /**
     * @brief Операция добавления сетки к существующей на основе элементарного объединения
     * @param mesh Указатель на сетку для добавления
     */
    Q_INVOKABLE void add(const TriangleMesh3D *mesh);
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
