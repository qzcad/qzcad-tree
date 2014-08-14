/**
  * Проект для консольного тестирования конечно элементных расчетов
  * */
#include <iostream>
#include <fstream>
#include <iomanip>
#include "hexahedralfem.h"

using namespace std;
using namespace msh;

/**
 * @brief isEquil - функция приблизительной проверки равенства двух действительных чисел
 * @param a - первое число
 * @param b - второе число
 * @param tol - точность
 * @return true, если a равно b, false - иначе
 */
inline bool isEquil (const Floating &a, const Floating &b, Floating tol = 0.00001)
{
    return fabs(a - b) < tol;
}

class XCondition: public FEMCondition3D
{
public:
    virtual bool isApplied(PointPointer point)
    {
        return isEquil(point->x(), 0.0);
    }
    virtual double u()
    {
        return 0.0;
    }
    virtual bool isU()
    {
        return true;
    }
    virtual double v()
    {
        return 0.0;
    }
    virtual bool isV()
    {
        return false;
    }
    virtual double w()
    {
        return 0.0;
    }
    virtual bool isW()
    {
        return false;
    }
};

class ZCondition: public FEMCondition3D
{
public:
    virtual bool isApplied(PointPointer point)
    {
        return isEquil(point->z(), 0.0);
    }
    virtual double u()
    {
        return 0.0;
    }
    virtual bool isU()
    {
        return false;
    }
    virtual double v()
    {
        return 0.0;
    }
    virtual bool isV()
    {
        return false;
    }
    virtual double w()
    {
        return 0.0;
    }
    virtual bool isW()
    {
        return true;
    }
};

class FixedCondition: public FEMCondition3D
{
public:
    virtual bool isApplied(PointPointer point)
    {
        return isEquil(point->x() * point->x() + point->z() * point->z(), 0.4 * 0.4);
    }
    virtual double u()
    {
        return 0.0;
    }
    virtual bool isU()
    {
        return true;
    }
    virtual double v()
    {
        return 0.0;
    }
    virtual bool isV()
    {
        return true;
    }
    virtual double w()
    {
        return 0.0;
    }
    virtual bool isW()
    {
        return true;
    }
};

class ForceCondition: public FEMCondition3D
{
public:
    virtual bool isApplied(PointPointer point)
    {
        return isEquil(point->y(), 0.018);
    }
    virtual double u()
    {
        return 0.0;
    }
    virtual bool isU()
    {
        return true;
    }
    virtual double v()
    {
        return -0.09;
    }
    virtual bool isV()
    {
        return true;
    }
    virtual double w()
    {
        return 0.0;
    }
    virtual bool isW()
    {
        return true;
    }
};

int main()
{
    cout << "Программа для тестирования расчетов при помощи МКЭ" << endl;
    fstream input;
    UInteger nodesCount;
    int freedom = 0;
    const int elementNodes = 8;
    HexahedralMesh3D mesh;
    UInteger elementsCount;

    const Floating G = 8.0E+4;
    const Floating nu = 0.27;
    const Floating E = 2.0 * G * (1.0 + nu);
    MechanicalParameters params(E, nu);

    cout << "Загрузка дискретной модели..."<< endl;
    input.open("plate_0_x_04_0_y_0018_0_z_04.txt", fstream::in);
    input >> freedom; // количество степеней сободы
    cout << "Степеней свободы - " << freedom << endl;
    input >> nodesCount;
    cout << "Узлов - " << nodesCount << "; загрузка узлов..." << endl;
    for (UInteger i = 0; i < nodesCount; i++)
    {
        Point3D node;
        Floating x, y, z;
        int type;
        input >> x;
        input >> y;
        input >> z;
        input >> type;
        node.set(x, y, z);
        mesh.pushNode(node, static_cast<NodeType>(type));
    }
    input >> elementsCount;
    cout << "Элементов - " << elementsCount << "; загрузка элементов..." << endl;
    for (UInteger i = 0; i < elementsCount; i++)
    {
        UInteger p[elementNodes];
        for (int j = 0; j < elementNodes; j++)
            input >> p[j];
        mesh.addElement(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]);
    }
    input.close();

    std::vector<FEMCondition3DPointer> boundary;
    XCondition xcond;
    ZCondition zcond;
    FixedCondition fixcond;
    boundary.push_back(&xcond);
    boundary.push_back(&zcond);
    boundary.push_back(&fixcond);
    ForceCondition force;

    HexahedralFEM fem(&mesh, params, &force, boundary);

    return 0;
}

