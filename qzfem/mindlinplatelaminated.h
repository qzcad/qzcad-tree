#ifndef MINDLINPLATELAMINATED_H
#define MINDLINPLATELAMINATED_H

#include "fem2d.h"
#include "femcondition.h"

#include "quadrilateralmesh2d.h"
#include "trianglemesh2d.h"
using namespace msh;

#include "doublematrix.h"
#include "doublevector.h"
using namespace mtx;


class MindlinPlateLaminated : public Fem2D
{
public:
    MindlinPlateLaminated(Mesh2D *mesh,
                          const std::vector<double> &thickness,
                          const std::vector<DoubleMatrix> &planeStressMatrix,
                          const std::list<FemCondition *> &conditions);
    MindlinPlateLaminated(Mesh2D *mesh,
                          const std::vector<double> &thickness,
                          const std::vector<DoubleMatrix> &D,
                          const std::vector<DoubleMatrix> &Dc,
                          const std::list<FemCondition *> &conditions);
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
protected:
    std::vector<double> thickness_; //!< Толщина объекта
    std::vector<DoubleMatrix> D_; //!< Матрица упругости
    std::vector<DoubleMatrix> Dc_; //!< Матрица упругости
};

#endif // MINDLINPLATELAMINATED_H
