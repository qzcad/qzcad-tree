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

class FixedCondition: public FEMCondition3D // защемление
{
public:
    virtual bool isApplied(PointPointer point)
    {
        return isEquil(point->x() * point->x() + point->z() * point->z(), 0.4 * 0.4);
//        return  isEquil(point->x() * point->x() + point->z() * point->z(), 0.2 * 0.2) || isEquil(point->x() * point->x() + point->z() * point->z(), 0.4 * 0.4);
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

class RotCondition: public FEMCondition3D // опирание
{
public:
    virtual bool isApplied(PointPointer point)
    {
        return  isEquil(point->y(), 0.0) && isEquil(point->x() * point->x() + point->z() * point->z(), 0.4 * 0.4);
//        return  isEquil(point->y(), 0.0) && isEquil(point->x() * point->x() + point->z() * point->z(), 0.2 * 0.2);
//        return  isEquil(point->y(), 0.0) && (isEquil(point->x() * point->x() + point->z() * point->z(), 0.2 * 0.2) || isEquil(point->x() * point->x() + point->z() * point->z(), 0.4 * 0.4));
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
        return true;
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
    int elementNodes = 8;
    HexahedralMesh3D mesh;
    UInteger elementsCount;
    int isNodeValue = 0;
    int isElementValue = 0;

    const Floating G = 8.0E+4;
    const Floating nu = 0.27;
    const Floating E = 2.0 * G * (1.0 + nu);
    MechanicalParameters3D params(E, nu);
    const Floating innerK = 6.0E+4;
    const Floating innerG = 2.77E+4;
    const Floating innerE = 9.0 * innerK * innerG / (3.0 * innerK + innerG);
    const Floating innerNu = (3.0 * innerK - 2.0 * innerG) / (2.0 * (3.0 * innerK + innerG));
    MechanicalParameters3D paramsInner(innerE, innerNu);

    cout << "Загрузка дискретной модели..."<< endl;
    input.open("plate_0_x_04_0_y_0018_0_z_04.txt", fstream::in);
    input >> freedom; // количество степеней сободы
    cout << "Степеней свободы - " << freedom << endl;
    input >> elementNodes;
    input >> nodesCount;
    cout << "Узлов - " << nodesCount << "; загрузка узлов..." << endl;
    input >> isNodeValue;
    for (UInteger i = 0; i < nodesCount; i++)
    {
        Point3D node;
        Floating x, y, z, val;
        int type;
        input >> x;
        input >> y;
        input >> z;
        input >> type;
        node.set(x, y, z);
        mesh.pushNode(node, static_cast<NodeType>(type));
        if (isNodeValue)
        {
            input >> val;
            mesh.pushNodeValue(val);
        }
    }
    input >> elementsCount;
    cout << "Элементов - " << elementsCount << "; загрузка элементов..." << endl;
    input >> isElementValue;
    for (UInteger i = 0; i < elementsCount; i++)
    {
        UInteger p[elementNodes];
        Floating val;
        for (int j = 0; j < elementNodes; j++)
            input >> p[j];
        mesh.addElement(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]);
        if (isElementValue)
        {
            input >> val;
            mesh.pushElementValue(val);
        }
    }
    input.close();

    std::vector<FEMCondition3DPointer> boundaryConditions;
    XCondition xcond;
    ZCondition zcond;
    FixedCondition fixcond;
    RotCondition rotcond;
    boundaryConditions.push_back(&xcond);
    boundaryConditions.push_back(&zcond);
//    boundaryConditions.push_back(&fixcond);
    boundaryConditions.push_back(&rotcond);
    ForceCondition force;
    std::vector<FEMCondition3DPointer> boundaryForces;
    boundaryForces.push_back(&force);

    HexahedralFEM fem(&mesh, params, boundaryForces, boundaryConditions);

    std::vector<MechanicalParameters3D> layers;
    layers.push_back(params);
    layers.push_back(paramsInner);
    layers.push_back(params);

//    HexahedralFEM fem(&mesh, layers, boundaryForces, boundaryConditions); // многослойный расчет

    fem.setNodeDisplacement(&mesh, 1);

    cout << "Сохранение результатов в файл" << endl;
    fstream output;
    output.open("v.txt", fstream::out);
    output << freedom << ' ' << elementNodes << endl;
    output << nodesCount << ' ' << 1 << endl;
    for (UInteger i = 0; i < nodesCount; i++)
    {
        PointPointer point = mesh.node(i);
        output << point->x() << ' ' << point->y() << ' ' << point->z() << ' ' << mesh.nodeType(i) << ' ' << mesh.nodeValue(i) << endl;
    }
    output << elementsCount << ' ' << 0 << endl;
    for (UInteger i = 0; i < elementsCount; i++)
    {
        ElementPointer element = mesh.element(i);
        for (int j = 0; j < element->verticesCount(); j++)
        {
            output << element->vertexNode(j) << " ";
        }
        output << endl;
    }
    output.close();

    return 0;
}

