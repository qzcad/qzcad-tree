#include "voxelmesh.h"

#include <iostream>

#include "consoleprogress.h"

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
    Point3D n(0.0, 0.0, 0.0);
    for (UInteger nnode: face)
    {
        n = n + normal_[nnode];
    }
    return n.normalized();
}

void VoxelMesh::voxel(const double &xMin, const double &yMin, const double &zMin, const double &width, const double &height, const double &depth, const double &h, std::function<double (double, double, double)> func, std::function<bool (double, double, double)> filter)
{
    clear();
    xMin_ = xMin; xMax_ = xMin + width;
    yMin_ = yMin; yMax_ = yMin + height;
    zMin_ = zMin; zMax_ = zMin + depth;
    if (xMin_ >= xMax_ || yMin_ >= yMax_ || zMin_ >= zMax_ || h <= 0.0)
    {
        std::cout << "Wrong dimensions: [" << xMin_ << "; " << xMax_ << "] x [" << yMin_ << "; " << yMax_ << "] x [" << zMin_ << "; " << zMax_ << "], h = " << h << std::endl;
        return;
    }
    ConsoleProgress progress(width / h);
    double x = xMin_;
    while (x < xMax_)
    {
        double y = yMin_;
        while (y < yMax_)
        {
            double z = zMin_;
            while (z < zMax_)
            {
                if (filter == nullptr || (filter != nullptr && filter(x + 0.5 * h, y + 0.5 * h, z + 0.5 * h)))
                {
                    addElement(addNode(Point3D(x, y, z), BORDER),
                               addNode(Point3D(x + h, y, z), BORDER),
                               addNode(Point3D(x + h, y + h, z), BORDER),
                               addNode(Point3D(x, y + h, z), BORDER),
                               addNode(Point3D(x, y, z + h), BORDER),
                               addNode(Point3D(x + h, y, z + h), BORDER),
                               addNode(Point3D(x + h, y + h, z + h), BORDER),
                               addNode(Point3D(x, y + h, z + h), BORDER));
                }
                z += h;
            }
            y += h;
        }
        x += h;
        ++progress;
    }
    progress.restart(nodesCount());
    for (Node3D n: node_)
    {
        Point3D g = grad(func, n.point, h);
        normal_.push_back(g);
        ++progress;
    }
    evalNodalValues(func);
}

}
