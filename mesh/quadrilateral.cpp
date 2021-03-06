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
    int d = i % 4;

    if (i < 0)
        return vertex_[(d == 0) ? d : d + 4];

    return vertex_[d];
}

const UInteger &Quadrilateral::operator [](int i) const
{
    int d = i % 4;

    if (i < 0)
        return vertex_[(d == 0) ? d : d + 4];

    return vertex_[d];
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
    int d = i % 4;

    if (i < 0)
        return vertex_[(d == 0) ? d : d + 4];

    return vertex_[d];
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

bool Quadrilateral::in(const UInteger &node) const
{
    return vertex_[0] == node || vertex_[1] == node || vertex_[2] == node || vertex_[3] == node;
}

int Quadrilateral::index(const UInteger &node) const
{
    if (vertex_[0] == node) return 0;
    if (vertex_[1] == node) return 1;
    if (vertex_[2] == node) return 2;
    if (vertex_[3] == node) return 3;
    return -1;
}
}
