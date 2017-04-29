#include "triangle.h"

namespace msh
{

Triangle::Triangle(const UInteger &node_0, const UInteger &node_1, const UInteger &node_2)
{
    vertex_[0] = node_0;
    vertex_[1] = node_1;
    vertex_[2] = node_2;
}

Triangle::Triangle(const Triangle &triangle)
{
    if (this != &triangle)
    {
        vertex_[0] = triangle.vertex_[0];
        vertex_[1] = triangle.vertex_[1];
        vertex_[2] = triangle.vertex_[2];
    }
}

UInteger &Triangle::operator [](int i)
{
    int d = i % 3;

    if (i < 0)
        return vertex_[(d == 0) ? d : d + 3];

    return vertex_[d];
}

const UInteger &Triangle::operator [](int i) const
{
    int d = i % 3;

    if (i < 0)
        return vertex_[(d == 0) ? d : d + 3];

    return vertex_[d];
}

Triangle &Triangle::operator =(const Triangle &triangle)
{
    if (this != &triangle)
    {
        vertex_[0] = triangle.vertex_[0];
        vertex_[1] = triangle.vertex_[1];
        vertex_[2] = triangle.vertex_[2];
    }
    return *this;
}

int Triangle::verticesCount() const
{
    return 3;
}

UInteger Triangle::vertexNode(int i) const
{
    int d = i % 3;

    if (i < 0)
        return vertex_[(d == 0) ? d : d + 3];

    return vertex_[d];
}

int Triangle::facesCount() const
{
    return 1;
}

UIntegerVector Triangle::face(const int &i) const
{
    UIntegerVector face(vertex_, vertex_ + 3);
    return face;
}

bool Triangle::in(const UInteger &node)
{
    return vertex_[0] == node || vertex_[1] == node || vertex_[2] == node;
}

int Triangle::index(const UInteger &node)
{
    if (vertex_[0] == node) return 0;
    if (vertex_[1] == node) return 1;
    if (vertex_[2] == node) return 2;
    return -1;
}

}
