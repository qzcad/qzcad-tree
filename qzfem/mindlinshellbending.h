#ifndef MINDLINSHELLBENDING_H
#define MINDLINSHELLBENDING_H

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
    MindlinShellBending(Mesh3D *mesh, double thickness, const DoubleMatrix &planeStressMatrix, const std::list<FemCondition *> &conditions, double alphaT = 0.0);
    MindlinShellBending(Mesh3D *mesh, double thickness, const DoubleMatrix &D, const DoubleMatrix &Dc, const std::list<FemCondition *> &conditions, double alphaT = 0.0);
    /**
     * @brief Конструктор для КЭ анаилза многослойных оболочек Миндлина
     * @param mesh Указатель на сетку элементов
     * @param thickness Массив толщин (для каждого слоя)
     * @param elasticMatrix Массив матриц упругих констант (для каждого слоя)
     * @param conditions Список условий (нагрузок и граничных условий)
     */
    //MindlinShellBending(Mesh3D *mesh, const std::vector<double> &thickness, const std::vector<DoubleMatrix> &planeStressMatrix, const std::list<FemCondition *> &conditions);
    /**
     * @brief Конструктор КЭ анализа разуршающих нагрузок для оболочки с использованием теории Миндлина и метода переменной жесткости
     * @param mesh Указатель на сетку элементов
     * @param thickness Толщина оболочки
     * @param strain Массив деформаций
     * @param stress Соответствующий массиву деформаций массив напряжений
     * @param nu Коэффициент Пуассона
     * @param conditions Список условий (нагрузок и граничных условий)
     */
   // MindlinShellBending(Mesh3D *mesh, double thickness, const std::vector<double> &strain, const std::vector<double> &stress, double nu, const std::list<FemCondition *> &conditions);
    MindlinShellBending(Mesh3D *mesh, std::function<double(double, double, double)> thickness_func, const DoubleMatrix &planeStressMatrix, const std::list<FemCondition *> &conditions, double alphaT = 0.0);
    MindlinShellBending(Mesh3D *mesh, std::function<double(double, double, double)> thickness_func, const DoubleMatrix &D, const DoubleMatrix &Dc, const std::list<FemCondition *> &conditions, double alphaT = 0.0);
    virtual void solve(std::function<double(double, double, double)> func, double delta=0.2, int maxiter=3);
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
    DoubleMatrix cosinuses(const Point3D &A, const Point3D &B, const Point3D &C);
    DoubleVector adaptationVector(const DoubleVector &displacement);
protected:
    double thickness_; //!< Толщина объекта
    DoubleMatrix D_; //!< Матрица упругости
    DoubleMatrix Dc_; //!< Матрица упругости
    double alpha_; //!< Коэффициент температурного напряжения
    std::function<double(double, double, double)> thickness_func_;
    DoubleVector mises_;
};

#endif // MINDLINSHELLBENDING_H
