#include "qdoublematrix.h"

QDoubleMatrix::QDoubleMatrix(QObject *parent) : QObject(parent)
{

}

QDoubleMatrix::QDoubleMatrix(const QDoubleMatrix &qdm) : QObject(qdm.parent()), DoubleMatrix(qdm)
{

}

QDoubleMatrix::QDoubleMatrix(const DoubleMatrix &dm): DoubleMatrix(dm)
{

}

QString QDoubleMatrix::toString() const
{
    QString s = "";
    for (size_type i = 0; i < rowCount(); i++)
    {
        s += " [ ";
        for (size_type j = 0; j < colCount(); j++)
        {
            s += tr("%1\t").arg(data(i, j));
        }
        s += "]\n";
    }
    return s;
}
