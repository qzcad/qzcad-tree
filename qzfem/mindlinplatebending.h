#ifndef MINDLINPLATEBENDING_H
#define MINDLINPLATEBENDING_H

#include <functional>

#include "fem2d.h"
#include "elasticmatrix.h"
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
    MindlinPlateBending(QuadrilateralMesh2D *mesh,
                        double thickness,
                        const ElasticMatrix &elasticMatrix,
                        std::list<FemCondition *> conditions);
    /**
     * @brief Метод возвращает название вектора узловых значений
     * @param num Номер вектора узловых значений
     * @return Название вектора узловых значений
     */
    virtual std::string nodeVectorName(UInteger num) const;
    /**
     * @brief Метод возвращает название вектора значений, определенных на элементе
     * @param num Номер вектора значений, определеных на элементе
     * @return Название вектора значений, определенных на элеменете
     */
    virtual std::string elementVectorName(UInteger num) const;
};

#endif // MINDLINPLATEBENDING_H
