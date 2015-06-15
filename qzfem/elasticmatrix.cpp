#include "elasticmatrix.h"

ElasticMatrix::ElasticMatrix(double E, double nu, bool isPlaneStress)
{
    DoubleMatrix d(3, 3);
    if (isPlaneStress)
    {
        d(0, 0) = 1.0;  d(0, 1) = nu;   d(0, 2) = 0.0;
        d(1, 0) = nu;   d(1, 1) = 1.0;  d(1, 2) = 0.0;
        d(2, 0) = 0.0;  d(2, 1) = 0.0;  d(2, 2) = (1 - nu) / 2.0;
        D_ = (E / (1.0 - nu*nu)) * d;
    }
    else
    {
        d(0, 0) = 1.0 - nu; d(0, 1) = nu;       d(0, 2) = 0.0;
        d(1, 0) = nu;       d(1, 1) = 1.0 - nu; d(1, 2) = 0.0;
        d(2, 0) = 0.0;      d(2, 1) = 0.0;      d(2, 2) = (1.0 - 2.0 * nu) / 2.0;
        D_ = (E / ((1.0 + nu) * (1.0 - 2.0 * nu))) * d;
    }
}

ElasticMatrix::ElasticMatrix(const ElasticMatrix &ematrix)
{
    D_ = ematrix.D_;
}
