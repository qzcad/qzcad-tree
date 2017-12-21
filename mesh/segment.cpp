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
    if ((i % 2) == 0)
        return leftVertex_;
    return rightVertex_;
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

UInteger &Segment::operator [](int i)
{
    if ((i % 2) == 0)
        return leftVertex_;
    return rightVertex_;
}

const UInteger &Segment::operator [](int i) const
{
    if ((i % 2) == 0)
        return leftVertex_;
    return rightVertex_;
}

bool Segment::isSame(const Segment &s)
{
    return (leftVertex_ == s.leftVertex_ && rightVertex_ == s.rightVertex_) ||
            (leftVertex_ == s.rightVertex_ && rightVertex_ == s.leftVertex_);
}

bool Segment::in(const UInteger &node) const
{
    return leftVertex_ == node || rightVertex_ == node;
}
}
