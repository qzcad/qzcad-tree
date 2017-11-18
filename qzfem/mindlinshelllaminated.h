#ifndef MINDLINSHELLLAMINATED_H
#define MINDLINSHELLLAMINATED_H
#include <funcopt.h>

#include "mindlinshellbending.h"

class MindlinShellLaminated : public MindlinShellBending
{
public:
    MindlinShellLaminated(Mesh3D *mesh, const std::vector<double> &thickness, const std::vector<DoubleMatrix> &planeStressMatrix, const std::list<FemCondition *> &conditions, double alphaT = 0.0);
    MindlinShellLaminated(Mesh3D *mesh, const std::vector<double> &thickness, const std::vector<DoubleMatrix> &D, const std::vector<DoubleMatrix> &Dc, const std::list<FemCondition *> &conditions, double alphaT = 0.0);
    MindlinShellLaminated(Mesh3D *mesh, std::function<std::vector<double>(double, double, double)> thickness, const std::vector<DoubleMatrix> &planeStressMatrix, const std::list<FemCondition *> &conditions, double alphaT = 0.0);
    MindlinShellLaminated(Mesh3D *mesh, std::function<std::vector<double>(double, double, double)> thickness, const std::vector<DoubleMatrix> &D, const std::vector<DoubleMatrix> &Dc, const std::list<FemCondition *> &conditions, double alphaT = 0.0);
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
    std::vector<double> thickness_; //!< Толщина объекта
    std::vector<DoubleMatrix> D_; //!< Матрица упругости
    std::vector<DoubleMatrix> Dc_; //!< Матрица упругости
    std::function<std::vector<double>(double, double, double)> thickness_func_;
};

#endif // MINDLINSHELLLAMINATED_H
