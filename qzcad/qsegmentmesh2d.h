#ifndef QSEGMENTMESH2D_H
#define QSEGMENTMESH2D_H

#include <QObject>
#include <qmetatype.h>
#include "segmentmesh2d.h"

using namespace msh;

/**
 * @brief Прокси-класс для создания сеток в скриптах моделей
 * @see QObject, QuadrilateralMesh2D
 */
class QSegmentMesh2D : public QObject, public SegmentMesh2D
{
    Q_OBJECT
public:
    /**
     * @brief Конструктор по умолчанию
     * @param parent Указатель на родительский объект
     */
    explicit QSegmentMesh2D(QObject *parent = 0);
    /**
     * @brief Конструктор копирования
     * @param qmesh Экземпляр объекта для копирования
     */
    QSegmentMesh2D(const QSegmentMesh2D &qmesh);
    /**
     * @brief Конструктор создает копию объекта, переданного по указателю
     * @param mesh Указатель на объект для копирования
     * @param parent Указатель на родительский объект
     */
    QSegmentMesh2D(SegmentMesh2D *mesh, QObject *parent = 0);
    /**
     * @brief Метод формирует строку с информацией о сетке
     * @return Строка с информацией о сетке
     * Метод доступен в скриптах
     */
    Q_INVOKABLE QString toString() const;
    Q_INVOKABLE double cfunction(const double &x, const double &y);
signals:

public slots:
};
// Регистрация класса для использования в скриптах моделей
Q_DECLARE_METATYPE(QSegmentMesh2D)
Q_DECLARE_METATYPE(QSegmentMesh2D*)
#endif // QSEGMENTMESH2D_H
