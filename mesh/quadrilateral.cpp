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
    if (i < 0)
        return vertex_[(i % 4) + 4];

    return vertex_[i % 4];
}

const UInteger &Quadrilateral::operator [](int i) const
{
    if (i < 0)
        return vertex_[(i % 4) + 4];

    return vertex_[i % 4];
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
    if (i < 0)
        return vertex_[(i % 4) + 4];

    return vertex_[i % 4];
}

int Quadrilateral::facesCount() const
{
    return 1;
}

UIntegerVector Quadrilateral::face(const int &i) const
{
    UIntegerVector face(vertex_, vertex_ + 4);
    return face;
}

bool Quadrilateral::in(const UInteger &node)
{
    return vertex_[0] == node || vertex_[1] == node || vertex_[2] == node || vertex_[3] == node;
}
}
