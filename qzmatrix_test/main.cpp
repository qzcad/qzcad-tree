#include <iostream>

using namespace std;

#include "doublevector.h"
#include "doublematrix.h"
#include "mappeddoublematrix.h"

using namespace mtx;

int main()
{
    DoubleVector dv(5, 0), dv1(10, 1);
    cout << dv.size() << endl;
    dv.print();
    dv1[5] = 5.5;
    dv1.print();
    dv = dv1;
    cout << dv.min() << "\t" << dv.max() << "\t" << dv.norm_2() << endl;
    dv.print();
    cout << "Hello World!" << endl;
    DoubleMatrix a(2, 0.0);
    a(0, 0) = 1.0;
    a(1, 1) = 2.5;
    a[0][1] = 0.3;
    a.print();
    DoubleMatrix b(2, 5, 1.0);
    b.print();
    cout << "a * b" << endl;
    DoubleMatrix m(a * b);
    m.print();
    b = a;
    b.print();
    DoubleMatrix c(b.transpose());
    c.print();
    cout << "sparce" << endl;
    MappedDoubleMatrix mdm(3);
    mdm(0, 1) = 0.5;
    mdm(1, 0) = 0.1; mdm(1, 2) = 0.3;
    mdm(2, 1) = 0.7;
    mdm.print();
    return 0;
}

