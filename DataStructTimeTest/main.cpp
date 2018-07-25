#include <iostream>

#include <time.h>

using namespace std;

#include "meshlist.h"
#include "hexahedralmesh3d.h"
using namespace msh;

#include "rfunctions.h"

double test(double x, double y, double z)
{
    return con(1.0 - x*x, con(1.0 - y*y, con(1.0 - z*z, con(x*x + y*y - 0.25, con(z*z + y*y - 0.25, con(z*z + x*x - 0.25, 1.3*1.3 - x*x - y*y - z*z))))));
}

int main()
{
    clock_t tStart;
    double elapsed_time;


    for (long d = 400; d <= 8000; d+=400)
    {
        int xCount = 30;
        int yCount = 30;
        int zCount = 30;
        double xMin = -1.5, yMin = -1.5, zMin = -1.5, width = 3.0, height = 3.0, depth = 3.0;
        double hx = width / (double)(xCount - 1);
        double hy = height / (double)(yCount - 1);
        double hz = depth / (double)(zCount - 1);
        float list_init=1.0e6, list_del=1.0e6, array_init=1.0e6, array_del=1.0e6;
        for (short k = 0; k < 3; k++)
        {
            MeshList ml;
//            cout << "Mesh List data structure" << endl;
//            cout << "Generating initial mesh... ";
            // формирование массива узлов
            tStart = clock();
            for (int i = 0; i < xCount - 1; i++)
            {
                double x = xMin + (double) i * hx;
                for (int j = 0; j < yCount - 1; j++)
                {
                    double y = yMin + (double) j * hy;
                    for (int k = 0; k < zCount - 1; k++)
                    {
                        double z = zMin + (double)k * hz;
                        Point3D p0(x, y, z);
                        Point3D p1(x + hx, y, z);
                        Point3D p2(x + hx, y + hy, z);
                        Point3D p3(x, y + hy, z);
                        Point3D p4(x, y, z + hz);
                        Point3D p5(x + hx, y, z + hz);
                        Point3D p6(x + hx, y + hy, z + hz);
                        Point3D p7(x, y + hy, z + hz);
                        std::vector<NodeList3d *> v (8);
                        v[0] = ml.addNode(p0, INNER);
                        v[1] = ml.addNode(p1, INNER);
                        v[2] = ml.addNode(p2, INNER);
                        v[3] = ml.addNode(p3, INNER);
                        v[4] = ml.addNode(p4, INNER);
                        v[5] = ml.addNode(p5, INNER);
                        v[6] = ml.addNode(p6, INNER);
                        v[7] = ml.addNode(p7, INNER);
                        ml.addElement(v);
                    }
                }
            }
            elapsed_time = (double)(clock() - tStart) / CLOCKS_PER_SEC;
            if (elapsed_time < list_init)
                list_init = elapsed_time;
//            cout <<"Done in " << elapsed_time << "; nodes: " << ml.nodes.size() << "; elements: " << ml.elements.size() << endl;
//            cout << "Deleting outside elements... ";
            tStart = clock();
            ml.clearFuncNodes(test, d);
            elapsed_time = (double)(clock() - tStart) / CLOCKS_PER_SEC;
            if (elapsed_time < list_del)
                list_del = elapsed_time;
//            cout <<"Done in " << elapsed_time << "; nodes: " << ml.nodes.size() << "; elements: " << ml.elements.size() << endl;

//            cout << "Mesh Array Data structure" << endl;
            HexahedralMesh3D hmesh;
//            cout << "Generating initial mesh... ";
            // формирование массива узлов
            tStart = clock();
            for (int i = 0; i < xCount - 1; i++)
            {
                double x = xMin + (double) i * hx;
                for (int j = 0; j < yCount - 1; j++)
                {
                    double y = yMin + (double) j * hy;
                    for (int k = 0; k < zCount - 1; k++)
                    {
                        double z = zMin + (double)k * hz;
                        Point3D p0(x, y, z);
                        Point3D p1(x + hx, y, z);
                        Point3D p2(x + hx, y + hy, z);
                        Point3D p3(x, y + hy, z);
                        Point3D p4(x, y, z + hz);
                        Point3D p5(x + hx, y, z + hz);
                        Point3D p6(x + hx, y + hy, z + hz);
                        Point3D p7(x, y + hy, z + hz);
                        std::vector<UInteger> v (8);
                        v[0] = hmesh.addNode(p0, INNER);
                        v[1] = hmesh.addNode(p1, INNER);
                        v[2] = hmesh.addNode(p2, INNER);
                        v[3] = hmesh.addNode(p3, INNER);
                        v[4] = hmesh.addNode(p4, INNER);
                        v[5] = hmesh.addNode(p5, INNER);
                        v[6] = hmesh.addNode(p6, INNER);
                        v[7] = hmesh.addNode(p7, INNER);
                        hmesh.addElement(v);
                    }
                }
            }
            elapsed_time = (double)(clock() - tStart) / CLOCKS_PER_SEC;
            if (elapsed_time < array_init)
                array_init = elapsed_time;
//            cout <<"Done in " << elapsed_time << "; nodes: " << hmesh.nodesCount() << "; elements: " << hmesh.elementsCount() << endl;
//            cout << "Deleting outside elements... ";
            tStart = clock();
            hmesh.clearFuncNodes(test, d);
            elapsed_time = (double)(clock() - tStart) / CLOCKS_PER_SEC;
            if (elapsed_time < array_del)
                array_del = elapsed_time;
//            cout <<"Done in " << elapsed_time << "; nodes: " << hmesh.nodesCount() << "; elements: " << hmesh.elementsCount() << endl;
        }
        cout << "<========================" << d << "============================>" << endl;
        cout << "List Mesh: init - " << list_init << " sec., del - " << list_del << " sec." << endl;
        cout << "Array Mesh: init - " << array_init << " sec., del - " << array_del << " sec." << endl;
    }
    cout << "Press any key...";
    return 0;
}
