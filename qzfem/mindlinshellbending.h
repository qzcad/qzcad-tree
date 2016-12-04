#ifndef MINDLINSHELLBENDING_H
#define MINDLINSHELLBENDING_H

#include <functional>
#include <list>

#include "fem2d.h"
#include "elasticmatrix.h"
#include "femcondition.h"

#include "quadrilateralmesh3d.h"
#include "trianglemesh3d.h"
using namespace msh;

#include "doublematrix.h"
#include "doublevector.h"
using namespace mtx;

class MindlinShellBending : public Fem2D
{
public:
    /**
     * @brief Конструктор для КЭ анализа однослойных оболочек Миндлина
     * @param mesh Указатель на сетку элементов
     * @param thickness Толщина оболочки
     * @param elasticMatrix Матрица упругих констант
     * @param conditions Список условий (нагрузок и граничных условий)
     */
    MindlinShellBending(Mesh3D *mesh, double thickness, const ElasticMatrix &elasticMatrix, std::list<FemCondition *> conditions);
    /**
     * @brief Конструктор для КЭ анаилза многослойных оболочек Миндлина
     * @param mesh Указатель на сетку элементов
     * @param thickness Массив толщин (для каждого слоя)
     * @param elasticMatrix Массив матриц упругих констант (для каждого слоя)
     * @param conditions Список условий (нагрузок и граничных условий)
     */
    MindlinShellBending(Mesh3D *mesh, const std::vector<double> &thickness, const std::vector<ElasticMatrix> &elasticMatrix, std::list<FemCondition *> conditions);
    /**
     * @brief Конструктор КЭ анализа разуршающих нагрузок для оболочки с использованием теории Миндлина и метода переменной жесткости
     * @param mesh Указатель на сетку элементов
     * @param thickness Толщина оболочки
     * @param strain Массив деформаций
     * @param stress Соответствующий массиву деформаций массив напряжений
     * @param nu Коэффициент Пуассона
     * @param conditions Список условий (нагрузок и граничных условий)
     */
    MindlinShellBending(Mesh3D *mesh, double thickness, const std::vector<double> &strain, const std::vector<double> &stress, double nu, std::list<FemCondition *> conditions);
private:
    DoubleMatrix cosinuses(const Point3D &A, const Point3D &B, const Point3D &C);
    DoubleVector evalForces(Mesh3D *mesh, std::list<FemCondition *> conditions);
    void processInitialValues(Mesh3D *mesh, std::list<FemCondition *> conditions, MappedDoubleMatrix &global, DoubleVector &force);
};

#endif // MINDLINSHELLBENDING_H
