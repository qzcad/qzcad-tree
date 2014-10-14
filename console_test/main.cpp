/**
  * Проект для консольного тестирования конечно элементных расчетов
  * */
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
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
inline bool isEquil (const double &a, const double &b, double tol = 0.00001)
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
        return -0.011;
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

class SpaceBoundaryZ: public FEMCondition3D
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

class SpaceBoundaryY: public FEMCondition3D
{
public:
    virtual bool isApplied(PointPointer point)
    {
        return isEquil(point->y(), 0.0);
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

class SpaceBoundaryX: public FEMCondition3D
{
public:
    virtual bool isApplied(PointPointer point)
    {
        return isEquil(point->x(), 4.014);
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

class SpaceForceX: public FEMCondition3D
{
public:
    virtual bool isApplied(PointPointer point)
    {
        return isEquil(point->x(), 0.0);
    }
    virtual double u()
    {
//        return -3.303567151;
        const double r = 1.99;
        const double l = 0.045;
        return 251.12 * 9806.65 * 1.0E-6 / (M_PI * (r*r - (r-l)*(r-l))); // MPa
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

void save_vector(const std::vector<double> &v, const char *name)
{
    fstream out;
    out.open(name, fstream::out);
    cout << name << ": " << v.size() << endl;
    out << v.size() << endl;
    for (UInteger i = 0; i < v.size(); i++)
    {
        out << std::fixed << std::setprecision(15) << v[i] << endl;
    }
    out.close();
}

int main()
{
    cout << "Программа для тестирования расчетов при помощи МКЭ" << endl;
    fstream input;
    UInteger nodesCount;
    int freedom = 0;
    int elementNodes = 8;
    HexahedralMesh3D mesh;
    UInteger elementsCount;
    int isLayers = 0;

    const double G = 8.0E+4;
    const double nu = 0.27;
    const double E = 2.0 * G * (1.0 + nu);
    MechanicalParameters3D params(E, nu);
    const double innerK = 6.0E+4;
    const double innerG = 2.77E+4;
    const double innerE = 9.0 * innerK * innerG / (3.0 * innerK + innerG);
    const double innerNu = (3.0 * innerK - 2.0 * innerG) / (2.0 * (3.0 * innerK + innerG));
    MechanicalParameters3D paramsInner(innerE, innerNu);

    cout << "Загрузка дискретной модели..."<< endl;
//    input.open("plate_0_x_04_0_y_0018_0_z_04.txt", fstream::in);
    input.open("spacecraft.txt", fstream::in);
    input >> freedom; // количество степеней сободы
    cout << "Степеней свободы - " << freedom << endl;
    input >> elementNodes;
    input >> nodesCount;
    cout << "Узлов - " << nodesCount << "; загрузка узлов..." << endl;
    for (UInteger i = 0; i < nodesCount; i++)
    {
        Point3D node;
        double x, y, z;
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
    input >> isLayers;
    mesh.clearLayers();
    for (UInteger i = 0; i < elementsCount; i++)
    {
        UInteger p[elementNodes];
        int val;
        for (int j = 0; j < elementNodes; j++)
            input >> p[j];
        mesh.addElement(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]);
        if (isLayers)
        {
            input >> val;
            mesh.pushLayer(val);
        }
    }
    input.close();

//    std::vector<FEMCondition3DPointer> boundaryConditions;
//    XCondition xcond;
//    ZCondition zcond;
//    FixedCondition fixcond;
//    RotCondition rotcond;
//    boundaryConditions.push_back(&xcond);
//    boundaryConditions.push_back(&zcond);
////    boundaryConditions.push_back(&fixcond);
//    boundaryConditions.push_back(&rotcond);
//    ForceCondition force;
//    std::vector<FEMCondition3DPointer> boundaryForces;
//    boundaryForces.push_back(&force);

////    HexahedralFEM fem(&mesh, params, boundaryForces, boundaryConditions);

//    std::vector<MechanicalParameters3D> layers;
//    layers.push_back(params);
//    layers.push_back(paramsInner);
//    layers.push_back(params);

//    HexahedralFEM fem(&mesh, layers, boundaryForces, boundaryConditions); // многослойный расчет


    // spacecraft
    std::vector<FEMCondition3DPointer> boundaryConditions;
    SpaceBoundaryX sbx;
    SpaceBoundaryY sby;
    SpaceBoundaryZ sbz;
    boundaryConditions.push_back(&sbx);
    boundaryConditions.push_back(&sby);
    boundaryConditions.push_back(&sbz);
    std::vector<FEMCondition3DPointer> forces;
    SpaceForceX sfx;
    forces.push_back(&sfx);
    std::vector<MechanicalParameters3D> layers;
    MechanicalParameters3D shpangout(50000.0, 0.3);
    MechanicalParameters3D penoplast(100.0, 0.2, 28.0);
    MechanicalParameters3D zapolnitel(15.0, 0.33, 350.0);
    MechanicalParameters3D obshivka(63000.0, 36000.0, 36000.0, 0.4, 0.4, 0.33, 63000.0 / (2.0 + 2.0 * 0.4), 36000.0 / (2.0 + 2.0 * 0.4), 36000.0 / (2.0 + 2.0 * 0.33));
    layers.push_back(shpangout); // 0
    layers.push_back(shpangout); // 1
    layers.push_back(shpangout); // 2
    layers.push_back(shpangout); // 3
    layers.push_back(penoplast); // 4
    layers.push_back(zapolnitel); // 5
    layers.push_back(obshivka); // 6
    layers.push_back(obshivka); // 7
    layers.push_back(obshivka); // 8
    layers.push_back(obshivka); // 9
    layers.push_back(zapolnitel); // 10
    layers.push_back(obshivka); // 11
    layers.push_back(shpangout); // 12
    layers.push_back(shpangout); // 13
    layers.push_back(shpangout); // 14
    layers.push_back(shpangout); // 15
    layers.push_back(penoplast); // 16
    layers.push_back(zapolnitel); // 17
    layers.push_back(obshivka); // 18
    layers.push_back(obshivka); // 19
    layers.push_back(obshivka); // 20

    HexahedralFEM fem(&mesh, layers, forces, boundaryConditions); // многослойный расчет


    cout << "Сохранение результатов в файл" << endl;
    save_vector(fem.u(), "u.txt");
    save_vector(fem.v(), "v.txt");
    save_vector(fem.w(), "w.txt");
    save_vector(fem.sigmaX(), "sigma_x.txt");
    save_vector(fem.sigmaY(), "sigma_y.txt");
    save_vector(fem.sigmaZ(), "sigma_z.txt");
    save_vector(fem.tauXY(), "tau_xy.txt");
    save_vector(fem.tauYZ(), "tau_yz.txt");
    save_vector(fem.tauZX(), "tau_zx.txt");
    save_vector(fem.sigma(), "sigma.txt");

    return 0;
}

