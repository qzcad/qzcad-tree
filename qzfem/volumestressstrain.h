#ifndef VOLUMESTRESSSTRAIN_H
#define VOLUMESTRESSSTRAIN_H

#include <functional>

#include "fem3d.h"

#include "mesh3d.h"
using namespace msh;

#include "doublematrix.h"
#include "doublevector.h"
using namespace mtx;

class VolumeStressStrain : public Fem3D
{
public:
    VolumeStressStrain(Mesh3D *mesh, const DoubleMatrix &D, const std::list<FemCondition *> &conditions);
    VolumeStressStrain(Mesh3D *mesh, std::function<DoubleMatrix(double, double, double)> func, const std::list<FemCondition *> &conditions);
    virtual void solve(std::function<double(double, double, double)> func, double delta=0.2, int maxiter=3);
    static DoubleMatrix evalStressMatrix(const double &E, const double &nu);
protected:
    /**
     * @brief Метод для построения глобальной матрицы системы
     */
    virtual void buildGlobalMatrix();
    /**
     * @brief Метод для построения вектора системы
     */
    virtual void buildGlobalVector();
    /**
     * @brief Метод для обработки результтов решения
     * @param nodalValues
     */
    virtual void processSolution(const DoubleVector &displacement);
    DoubleVector adaptationVector(const DoubleVector &displacement);
protected:
    DoubleMatrix D_; //!< Матрица упругости
    std::function<DoubleMatrix(double, double, double)> D_func_; //!< функциональная матрица упругости
};

#endif // VOLUMESTRESSSTRAIN_H
