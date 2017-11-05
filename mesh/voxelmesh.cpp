#include "voxelmesh.h"

#include <iostream>

namespace msh {

VoxelMesh::VoxelMesh(): HexahedralMesh3D()
{

}

VoxelMesh::VoxelMesh(const VoxelMesh &mesh): HexahedralMesh3D(mesh)
{
    normal_ = mesh.normal_;
}

VoxelMesh::VoxelMesh(const VoxelMesh *mesh): HexahedralMesh3D(mesh)
{
    normal_ = mesh->normal_;
}

void VoxelMesh::clearElements()
{
    normal_.clear();
    HexahedralMesh3D::clearElements();
}

Point3D VoxelMesh::normal(const UIntegerVector &face) const
{
    Point3D n = normal_[face[0]];
    for (int i = 1; i < face.size(); i++)
    {
        n = n + normal_[face[i]];
    }
    return n.normalized();
}

void VoxelMesh::voxel(const double &xMin, const double &yMin, const double &zMin, const double &width, const double &height, const double &depth, const double &h, std::function<double (double, double, double)> func, std::function<bool (double, double, double)> filter)
{
    clear();
    xMin_ = xMin; xMax_ = xMin + width;
    yMin_ = yMin; yMax_ = yMin + height;
    zMin_ = zMin; zMax_ = zMin + depth;
    if (xMin_ >= xMax_ || yMin_ >= yMax_ || zMin_ >= zMax_ || h < 0)
    {
        std::cout << "Wrong dimensions: [" << xMin_ << "; " << xMax_ << "] x [" << yMin_ << "; " << yMax_ << "] x [" << zMin_ << "; " << zMax_ << "], h = " << h << std::endl;
        return;
    }
    double x = xMin_;
    while (x < xMax_)
    {
        double y = yMin_;
        while (y < yMax_)
        {
            double z = zMin_;
            while (z < zMax_)
            {
                if (filter == nullptr ||
                        (filter != nullptr && filter(x, y, z) && filter(x + h, y, z) && filter(x + h, y + h, z) && filter(x, y + h, z)
                         && filter(x, y, z + h) && filter(x + h, y, z + h) && filter(x + h, y + h, z + h) && filter(x, y + h, z + h)))
                {
                    addElement(pushNode(Point3D(x, y, z), BORDER),
                               pushNode(Point3D(x + h, y, z), BORDER),
                               pushNode(Point3D(x + h, y + h, z), BORDER),
                               pushNode(Point3D(x, y + h, z), BORDER),
                               pushNode(Point3D(x, y, z + h), BORDER),
                               pushNode(Point3D(x + h, y, z + h), BORDER),
                               pushNode(Point3D(x + h, y + h, z + h), BORDER),
                               pushNode(Point3D(x, y + h, z + h), BORDER));
                }
                z += h;
            }
            y += h;
        }
        x += h;
    }
    for (Node3D n: node_)
    {
        Point3D g = grad(func, n.point, h);
        normal_.push_back(g);
    }
    evalNodalValues(func);
}

}
