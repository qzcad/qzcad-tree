#ifndef MINDLINSHELLLAMINATED_H
#define MINDLINSHELLLAMINATED_H

#include "mindlinshellbending.h"

class MindlinShellLaminated : public MindlinShellBending
{
public:
    MindlinShellLaminated(Mesh3D *mesh, const std::vector<double> &thickness, const std::vector<DoubleMatrix> &planeStressMatrix, const std::list<FemCondition *> &conditions, double alphaT = 0.0);
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
};

#endif // MINDLINSHELLLAMINATED_H
