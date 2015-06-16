/**
  * @author Сергей Чопоров
  * @date 15/06/2015
  * @version 1.0.5
  * @copyright Copyright 2015 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  * */
#ifndef QUADRILATERALFEM_H
#define QUADRILATERALFEM_H

#include <functional>
#include <vector>

#include "quadrilateralmesh2d.h"
#include "elasticmatrix.h"

using namespace msh;

/**
 * @brief Класс для решения задач МКЭ на базе четырехугольных конечных элементов
 */
class QuadrilateralFEM
{
public:
    /**
     * @brief Конструктор по умолчанию
     */
    QuadrilateralFEM();
    void planeStressStrain(QuadrilateralMesh2D *mesh,
                           double thickness,
                           const ElasticMatrix &elasticMatrix,
                           std::function<int(double, double)> fixFunc,
                           std::function<Point2D(double, double)> boundaryValue,
                           std::function<Point2D(double, double)> nodalForce,
                           std::function<Point2D(double, double)> surfaceForce,
                           std::function<Point2D(double, double)> volumeForce);
    std::vector<double> u() const;
    void setU(const std::vector<double> &u);

    std::vector<double> v() const;
    void setV(const std::vector<double> &v);

    std::vector<double> sigmaX() const;
    void setSigmaX(const std::vector<double> &sigmaX);

    std::vector<double> sigmaY() const;
    void setSigmaY(const std::vector<double> &sigmaY);

    std::vector<double> tauXY() const;
    void setTauXY(const std::vector<double> &tauXY);

private:
    std::vector<double> u_; //!< Перемещения в первом направлении (x)
    std::vector<double> v_; //!< Перемещения во втором направлении (y)
    std::vector<double> sigmaX_; //!< Нормальные компоненты напряжения, параллельные первому направлению (x)
    std::vector<double> sigmaY_; //!< Нормальные компоненты напряжения, параллельные второму направлению (y)
    std::vector<double> tauXY_; //!< Касательные компоненты напряжения, в плоскости 1-2 (x-y)
};

#endif // QUADRILATERALFEM_H
