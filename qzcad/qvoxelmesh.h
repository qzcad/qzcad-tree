#ifndef QVOXELMESH_H
#define QVOXELMESH_H

#include <QObject>
#include <qmetatype.h>

#include "voxelmesh.h"

using namespace msh;

class QVoxelMesh: public QObject, public VoxelMesh
{
    Q_OBJECT
public:
    /**
     * @brief Конструктор по умолчанию
     * @param parent Указатель на родительский объект
     */
    explicit QVoxelMesh(QObject *parent = 0);
    /**
     * @brief Конструктор копирования
     * @param qmesh Экземпляр объекта для копирования
     */
    QVoxelMesh(const QVoxelMesh &qmesh);
    /**
     * @brief Конструктор создает копию объекта, переданного по указателю
     * @param mesh Указатель на объект для копирования
     * @param parent Указатель на родительский объект
     */
    QVoxelMesh(QVoxelMesh *mesh, QObject *parent = 0);
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
Q_DECLARE_METATYPE(QVoxelMesh)
Q_DECLARE_METATYPE(QVoxelMesh*)

#endif // QVOXELMESH_H
