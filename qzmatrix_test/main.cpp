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
    cout << "a:" << endl;
    a.print();
    cout << "|a| = " << a.det2x2() << endl;
    DoubleMatrix e = a * a.inverted2x2();
    cout << "e:" << endl;
    e.print();
    DoubleMatrix b(2, 5, 1.0);
    cout << "b:" << endl;
    b.print();

    DoubleMatrix m(a * b);
    cout << "a * b" << endl;
    m.print();
    b = a;
    cout << "b = a :" << endl;
    b.print();
    DoubleMatrix c(b.transpose());
    cout << "c = b.traspose():" << endl;
    c.print();
    DoubleMatrix ccc(3, 0.0);
    ccc(0, 1) = 0.8;
    ccc(1, 0) = 0.2; ccc(1, 2) = 0.4;
    ccc(2, 1) = 0.5;
    ccc(0, 2) = 1.0;
    cout << "ccc:" << endl;
    ccc.print();
    cout << "|ccc| = " << ccc.det3x3() << endl;
    DoubleMatrix iccc = ccc.inverted3x3();
    cout << "ccc^(-1):" << endl;
    iccc.print();
    DoubleMatrix eee = ccc * iccc;
    cout << "eee = ccc * ccc^(-1) :" << endl;
    eee.print();
    cout << "sparce" << endl;
    MappedDoubleMatrix mdm(3);
//    mdm(0, 1) = 0.5;
//    mdm(1, 0) = 0.1; mdm(1, 2) = 0.3;
//    mdm(2, 1) = 0.7;
    mdm(0, 0) = 1.0; mdm(0, 1) = 1.0;
    mdm(1, 0) = 1.0; mdm(1, 1) = 1.0; mdm(1, 2) = 2.0;
                     mdm(2, 1) = 2.0; mdm(2, 2) = 1.0;
    mdm.print();
    DoubleVector v3(3, 1.0);
    v3[0] += 1.0;
    cout << "v3: ";
    v3.print();
    DoubleVector ccc_v3(mdm * v3);
    cout << "mdm * v3: ";
    ccc_v3.print();
    cout << "v * ccc_v3 = " << v3 * ccc_v3 << endl;
    DoubleVector X = mdm.conjugateGradient(v3);
    X.print();
//    MappedDoubleMatrix tst(2);
//    tst(0,0) = 4.0; tst(0, 1) = 1.0;
//    tst(1,0) = 1.0; tst(1, 1) = 3.0;
//    cout << "test system: " << endl;
//    tst.print();
//    DoubleVector tsv(2);
//    tsv[0] = 1.0; tsv[1] = 2.0;
//    cout << endl;
//    tsv.print();
//    DoubleVector X = tst.conjugateGradient(tsv);
//    X.print();
    return 0;
}

