#include "quadrilateral.h"

namespace msh
{
Quadrilateral::Quadrilateral(const UInteger &node_0, const UInteger &node_1, const UInteger &node_2, const UInteger &node_3)
{
    vertex_[0] = node_0;
    vertex_[1] = node_1;
    vertex_[2] = node_2;
    vertex_[3] = node_3;
}

Quadrilateral::Quadrilateral(const Quadrilateral &quad)
{
    if (this != &quad)
    {
        vertex_[0] = quad.vertex_[0];
        vertex_[1] = quad.vertex_[1];
        vertex_[2] = quad.vertex_[2];
        vertex_[3] = quad.vertex_[3];
    }
}

UInteger &Quadrilateral::operator [](int i)
{
    switch (i)
    {
    case -1:
        return vertex_[3];
    case 0:
    case 1:
    case 2:
    case 3:
        return vertex_[i];
    case 4:
        return vertex_[0];
    }
    return vertex_[0];
}

const UInteger &Quadrilateral::operator [](int i) const
{
    switch (i)
    {
    case -1:
        return vertex_[3];
    case 0:
    case 1:
    case 2:
    case 3:
        return vertex_[i];
    case 4:
        return vertex_[0];
    }
    return vertex_[0];
}

Quadrilateral &Quadrilateral::operator =(const Quadrilateral &quad)
{
    if (this != &quad)
    {
        vertex_[0] = quad.vertex_[0];
        vertex_[1] = quad.vertex_[1];
        vertex_[2] = quad.vertex_[2];
        vertex_[3] = quad.vertex_[3];
    }
    return *this;
}

int Quadrilateral::verticesCount() const
{
    return 4;
}

UInteger Quadrilateral::vertexNode(int i) const
{
    switch (i)
    {
    case -1:
        return vertex_[3];
    case 0:
    case 1:
    case 2:
    case 3:
        return vertex_[i];
    case 4:
        return vertex_[0];
    }
    return vertex_[0];
}

int Quadrilateral::facesCount() const
{
    return 1;
}

UIntegerVector Quadrilateral::face(const int &i) const
{
    UIntegerVector face(vertex_, vertex_ + sizeof(vertex_) / sizeof(UInteger));
    return face;
}
}
