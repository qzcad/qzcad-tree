#include "tetrahedron.h"

namespace msh {

Tetrahedron::Tetrahedron(const UInteger &node_0, const UInteger &node_1, const UInteger &node_2, const UInteger &node_3)
{
    vertex_[0] = node_0;
    vertex_[1] = node_1;
    vertex_[2] = node_2;
    vertex_[3] = node_3;
}

Tetrahedron::Tetrahedron(const Tetrahedron &tetrahedron)
{
    vertex_[0] = tetrahedron.vertex_[0];
    vertex_[1] = tetrahedron.vertex_[1];
    vertex_[2] = tetrahedron.vertex_[2];
    vertex_[3] = tetrahedron.vertex_[3];
}

UInteger &Tetrahedron::operator [](int i)
{
    int d = i % 4;

    if (i < 0)
        return vertex_[(d == 0) ? d : d + 4];

    return vertex_[d];
}

Tetrahedron &Tetrahedron::operator =(const Tetrahedron &tetrahedron)
{
    if (this != &tetrahedron)
    {
        vertex_[0] = tetrahedron.vertex_[0];
        vertex_[1] = tetrahedron.vertex_[1];
        vertex_[2] = tetrahedron.vertex_[2];
        vertex_[3] = tetrahedron.vertex_[3];
    }
    return *this;
}

int Tetrahedron::verticesCount() const
{
    return 4;
}

UInteger Tetrahedron::vertexNode(int i) const
{
    int d = i % 4;

    if (i < 0)
        return vertex_[(d == 0) ? d : d + 4];

    return vertex_[d];
}

int Tetrahedron::facesCount() const
{
    return 4;
}

UIntegerVector Tetrahedron::face(const int &i) const
{
    UIntegerVector face(3);
    if (i == 0)
    {
        face[0] = vertex_[0];
        face[1] = vertex_[2];
        face[2] = vertex_[1];
    }
    else if (i == 1)
    {
        face[0] = vertex_[0];
        face[1] = vertex_[1];
        face[2] = vertex_[3];
    }
    else if (i == 2)
    {
        face[0] = vertex_[1];
        face[1] = vertex_[2];
        face[2] = vertex_[3];
    }
    else
    {
        face[0] = vertex_[2];
        face[1] = vertex_[0];
        face[2] = vertex_[3];
    }
    return face;
}

bool Tetrahedron::in(const UInteger &node) const
{
    return vertex_[0] == node || vertex_[1] == node || vertex_[2] == node || vertex_[3] == node;
}

}
