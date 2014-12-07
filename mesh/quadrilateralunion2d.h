/**
  * @author Сергей Чопоров
  * @date 19/03/2014
  * @version 1.0.0
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef QUADRILATERALUNION2D_H
#define QUADRILATERALUNION2D_H

#include "quadrilateralmesh2d.h"

namespace msh
{
/**
 * @brief Клас QuadrilateralUnion2D является абстракцией "простого" объединения сеток четырехугольных элементов
 * Под "простым" объединением понимается объединение без проверки наложения элементов.
 * Совпадающие узлы объединяются.
 * Для формирования массива признаков граничных узлов необходимо выполнить дополнительную процедуру.
 */
class QuadrilateralUnion2D : public QuadrilateralMesh2D
{
public:
    QuadrilateralUnion2D();
    void addQuadRegion(const UInteger &xCount, const UInteger &yCount, const Point2D &v0, const Point2D &v1, const Point2D &v2, const Point2D &v3);
    void addTriangleRegion(const UInteger &count, const Point2D &v0, const Point2D &v1, const Point2D &v2);
    void addMesh(QuadrilateralMesh2D *mesh);
private:
    bool isFirstRegion_; //!< Признак добавления первой элементарной области
    int layerNumber_;
};
}

#endif // QUADRILATERALUNION2D_H
