#ifndef QDOUBLEMATRIX_H
#define QDOUBLEMATRIX_H

#include <QObject>
#include <qmetatype.h>
#include "doublematrix.h"

using namespace mtx;

class QDoubleMatrix : public QObject, public DoubleMatrix
{
    Q_OBJECT
public:
    explicit QDoubleMatrix(QObject *parent = 0);
    QDoubleMatrix(const QDoubleMatrix &qdm);
    QDoubleMatrix(const DoubleMatrix &dm);
    Q_INVOKABLE QString toString() const;
};
Q_DECLARE_METATYPE(QDoubleMatrix)
Q_DECLARE_METATYPE(QDoubleMatrix*)
#endif // QDOUBLEMATRIX_H
