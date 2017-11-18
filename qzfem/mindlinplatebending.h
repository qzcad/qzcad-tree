#ifndef MINDLINPLATEBENDING_H
#define MINDLINPLATEBENDING_H

#include <functional>

#include "fem2d.h"
#include "femcondition.h"

#include "quadrilateralmesh2d.h"
#include "trianglemesh2d.h"
using namespace msh;

#include "doublematrix.h"
#include "doublevector.h"
using namespace mtx;

class MindlinPlateBending : public Fem2D
{
public:
    MindlinPlateBending(Mesh2D *mesh,
                        double thickness,
                        const DoubleMatrix &D,
                        const std::list<FemCondition *> &conditions);
    MindlinPlateBending(Mesh2D *mesh,
                        double thickness,
                        const DoubleMatrix &D,
                        const DoubleMatrix &Dc,
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
    double thickness_; //!< Толщина объекта
    DoubleMatrix D_; //!< Матрица упругости
    DoubleMatrix Dc_; //!< Матрица упругости для сдвиговых деформаций
};

#endif // MINDLINPLATEBENDING_H
