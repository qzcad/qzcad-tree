#ifndef MINDLINSHELLPLASTIC_H
#define MINDLINSHELLPLASTIC_H

#include <functional>
#include <list>

#include "fem2d.h"
#include "femcondition.h"

#include "quadrilateralmesh3d.h"
#include "trianglemesh3d.h"
using namespace msh;

#include "doublematrix.h"
#include "doublevector.h"
using namespace mtx;

#include "mindlinshellbending.h"

class MindlinShellPlastic : public MindlinShellBending
{
public:
    MindlinShellPlastic(Mesh3D *mesh, double thickness, std::function<double(double)> E, const double &sigma_yield, const double  &sigma_max, const double &nu, const std::list<FemCondition *> &conditions, double alphaT = 0.0);
    MindlinShellPlastic(Mesh3D *mesh, std::function<double(double, double, double)> thickness_func, std::function<double(double)> E, const double &sigma_yield, const double  &sigma_max, const double &nu, const std::list<FemCondition *> &conditions, double alphaT = 0.0);
    virtual void solve();
protected:
    /**
     * @brief Метод для построения глобальной матрицы системы
     */
    virtual void buildGlobalMatrix();
    /**
     * @brief Метод для обработки результтов решения
     * @param nodalValues
     */
    virtual void processSolution(const DoubleVector &displacement);
protected:
    std::function<double(double)> E_;
    double nu_;
    double sigma_yield_;
    DoubleVector x_;
    DoubleVector y_;
    DoubleVector z_;
    DoubleVector theta_x_;
    DoubleVector theta_y_;
    DoubleVector theta_z_;
    DoubleVector sigma_x_;
    DoubleVector sigma_y_;
    DoubleVector sigma_z_;
    DoubleVector tau_xy_;
    DoubleVector tau_xz_;
    DoubleVector tau_yz_;
    DoubleVector mises_;
    double sigma_max_;
    double f_;
};

#endif // MINDLINSHELLPLASTIC_H
