#include "hexahedral.h"

namespace msh {

Hexahedral::Hexahedral(const UInteger &node_0, const UInteger &node_1, const UInteger &node_2, const UInteger &node_3, const UInteger &node_4, const UInteger &node_5, const UInteger &node_6, const UInteger &node_7)
{
    vertex_[0] = node_0;
    vertex_[1] = node_1;
    vertex_[2] = node_2;
    vertex_[3] = node_3;
    vertex_[4] = node_4;
    vertex_[5] = node_5;
    vertex_[6] = node_6;
    vertex_[7] = node_7;
}

Hexahedral::Hexahedral(const Hexahedral &hexahedral)
{
    vertex_[0] = hexahedral.vertex_[0];
    vertex_[1] = hexahedral.vertex_[1];
    vertex_[2] = hexahedral.vertex_[2];
    vertex_[3] = hexahedral.vertex_[3];
    vertex_[4] = hexahedral.vertex_[4];
    vertex_[5] = hexahedral.vertex_[5];
    vertex_[6] = hexahedral.vertex_[6];
    vertex_[7] = hexahedral.vertex_[7];
}

UInteger &Hexahedral::operator [](int i)
{
    int d = i % 8;

    if (i < 0)
        return vertex_[(d == 0) ? d : d + 8];

    return vertex_[d];
}

Hexahedral &Hexahedral::operator =(const Hexahedral &hexahedral)
{
    if (this != &hexahedral)
    {
        vertex_[0] = hexahedral.vertex_[0];
        vertex_[1] = hexahedral.vertex_[1];
        vertex_[2] = hexahedral.vertex_[2];
        vertex_[3] = hexahedral.vertex_[3];
        vertex_[4] = hexahedral.vertex_[4];
        vertex_[5] = hexahedral.vertex_[5];
        vertex_[6] = hexahedral.vertex_[6];
        vertex_[7] = hexahedral.vertex_[7];
    }
    return *this;
}

int Hexahedral::verticesCount() const
{
    return 8;
}

UInteger Hexahedral::vertexNode(int i) const
{
    int d = i % 8;

    if (i < 0)
        return vertex_[(d == 0) ? d : d + 8];

    return vertex_[d];
}

int Hexahedral::facesCount() const
{
    return 6;
}

UIntegerVector Hexahedral::face(const int &i) const
{
    UIntegerVector face(4);
    if (i == 0)
    {   // 0_1_2_3
        face[0] = vertex_[0];   face[3] = vertex_[1];   face[2] = vertex_[2]; face[1] = vertex_[3];
    }
    else if (i == 1)
    {   // 0_4_5_1
        face[0] = vertex_[0];   face[3] = vertex_[4];   face[2] = vertex_[5]; face[1] = vertex_[1];
    }
    else if (i == 2)
    {   // 1_5_6_2
        face[0] = vertex_[1];   face[3] = vertex_[5];   face[2] = vertex_[6]; face[1] = vertex_[2];
    }
    else if (i == 3)
    {   // 3_2_6_7
        face[0] = vertex_[3];   face[3] = vertex_[2];   face[2] = vertex_[6]; face[1] = vertex_[7];
    }
    else if (i == 4)
    {   // 0_3_7_4
        face[0] = vertex_[0];   face[3] = vertex_[3];   face[2] = vertex_[7]; face[1] = vertex_[4];
    }
    else if (i == 5)
    {   // 4_7_6_5
        face[0] = vertex_[4];   face[3] = vertex_[7];   face[2] = vertex_[6]; face[1] = vertex_[5];
    }
    return face;
}

bool Hexahedral::in(const UInteger &node) const
{
    return vertex_[0] == node || vertex_[1] == node || vertex_[2] == node || vertex_[3] == node ||
            vertex_[4] == node || vertex_[5] == node || vertex_[6] == node || vertex_[7] == node;
}

}
