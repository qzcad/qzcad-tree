#ifndef VOXELMESH_H
#define VOXELMESH_H

#include "hexahedralmesh3d.h"

namespace msh {

class VoxelMesh : public HexahedralMesh3D
{
public:
    VoxelMesh();
    VoxelMesh(const VoxelMesh &mesh);
    VoxelMesh(const VoxelMesh *mesh);
    /**
     * @brief Метод очищает информацию об елементах
     */
    virtual void clearElements();
    /**
     * @brief Вычислить нормаль к грани
     * @param face Масив, представляющий грань
     * @return Координаты вектора-нормали
     */
    virtual Point3D normal(const UIntegerVector &face) const;
    void voxel(const double &xMin, const double &yMin, const double &zMin,
               const double &width, const double &height, const double &depth, const double &h,
               std::function<double(double, double, double)> func,
               std::function<bool(double, double, double)> filter = nullptr);
protected:
    std::vector<Point3D> normal_;
};

}

#endif // VOXELMESH_H
