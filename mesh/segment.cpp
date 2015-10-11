#include "segment.h"

namespace msh
{
Segment::Segment(const long &left, const long &right)
{
    leftVertex_ = left;
    rightVertex_ = right;
}

Segment::Segment(const Segment &segment)
{
    leftVertex_ = segment.leftVertex_;
    rightVertex_ = segment.rightVertex_;
}

int Segment::verticesCount() const
{
    return 2;
}

UInteger Segment::vertexNode(int i) const
{
    if (i == 0)
        return leftVertex_;
    if (i == 1)
        return rightVertex_;
    return leftVertex_;
}

int Segment::facesCount() const
{
    return 1;
}

UIntegerVector Segment::face(const int &i) const
{
    UIntegerVector f(2);
    f[0] = leftVertex_;
    f[1] = rightVertex_;
    return f;
}
}
