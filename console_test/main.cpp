/**
  * Проект для консольного тестирования конечно элементных расчетов
  * */
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "hexahedralfem.h"
#include "plasticfem.h"

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

class ForceCondition: public ForceCondition3D
{
public:
    ForceCondition()
    {
        setForceType(SURFACE_FORCE);
    }

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
        return -force_; // давление напралено в низ
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
    double force() const
    {
    return force_;
    }
    void setForce(double force)
    {
        force_ = force;
    }
private:
    double force_;
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

class SpaceForceX: public ForceCondition3D
{
public:
    SpaceForceX()
    {
        setForceType(SURFACE_FORCE);
    }

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

class TankBoundaryFixed: public FEMCondition3D
{
public:
    virtual bool isApplied(PointPointer point)
    {
        return isEquil(point->x(), 0) || isEquil(point->x(), 2.949);
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

class TankBoundaryXMove: public FEMCondition3D
{
public:
    virtual bool isApplied(PointPointer point)
    {
        return isEquil(point->x(), 0) || isEquil(point->x(), 2.949);
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
        return true;
    }
};

class TankForce: public ForceCondition3D
{
public:
    TankForce()
    {
        setForceType(SURFACE_FORCE);
    }

    virtual bool isApplied(PointPointer point)
    {
        x = point->x();
        y = point->y();
        z = point->z();
        return isEquil(y * y + z * z, 1.5 * 1.5);
    }
    virtual double u()
    {
        return 0.0; // MPa
    }
    virtual bool isU()
    {
        return true;
    }
    virtual double v()
    {
        return pressure * sin(atan2(y,z));
    }
    virtual bool isV()
    {
        return false;
    }
    virtual double w()
    {
        return pressure * cos(atan2(y,z));
    }
    virtual bool isW()
    {
        return false;
    }
    void setPressure(double p)
    {
        pressure = 0.101325 * p; // атмосферы в мегапаскали
    }

private:
    double x;
    double y;
    double z;
    double pressure;
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

    HexahedralFEM *fem = 0;

    int task = 0;

    cout << "Введите номер тестовой задачи:" << endl <<
            "1 - пластинка (однослойный случай);" << endl <<
            "2 - пластинка (три слоя);" << endl <<
            "3 - многослойная оболочка 2Ц41;" << endl <<
            "4 - цилиндрическая оболочка бака (упругий случай);" << endl <<
            "5 - цилиндрическая оболочка (упруго-пластичность)." << endl;
    cin >> task;

    if (task < 1 || task > 5)
    {
        cout << "Ошибка: теста с указанным номером не сучествует.";
        return 0;
    }

    cout << "Загрузка дискретной модели..."<< endl;

    if (task == 1 || task == 2)
        input.open("plate_0_x_04_0_y_0018_0_z_04.txt", fstream::in);
    else if (task == 3)
        input.open("spacecraft.txt", fstream::in);
    else if (task == 4 || task == 5)
        input.open("tank.txt", fstream::in);

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

    if (task == 1 || task == 2)
    {
        const double G = 8.0E+4;
        const double nu = 0.27;
        const double E = 2.0 * G * (1.0 + nu);
        MechanicalParameters3D params(E, nu);
        const double innerK = 6.0E+4;
        const double innerG = 2.77E+4;
        const double innerE = 9.0 * innerK * innerG / (3.0 * innerK + innerG);
        const double innerNu = (3.0 * innerK - 2.0 * innerG) / (2.0 * (3.0 * innerK + innerG));
        MechanicalParameters3D paramsInner(innerE, innerNu);
        std::vector<FEMCondition3DPointer> boundaryConditions;
        XCondition xcond;
        ZCondition zcond;
        FixedCondition fixcond;
        RotCondition rotcond;
        int boundary_type = 0;
        cout << "Выберите тип граничных условий (1 - защемеление; 2 - свободное обперание): ";
        cin >> boundary_type;
        if (boundary_type != 1 && boundary_type != 2)
        {
            cout << "Ошибка: некорректный код типа граничных условий." << endl;
            return 0;
        }
        boundaryConditions.push_back(&xcond);
        boundaryConditions.push_back(&zcond);
        if (boundary_type == 1)
            boundaryConditions.push_back(&fixcond);
        else
            boundaryConditions.push_back(&rotcond);
        ForceCondition force;
        double p = 0;
        cout << "Введите значение распределенного по верхней грани пластике давления (p > 0): ";
        cin >> p;
        if (p <= 0)
        {
            cout << "Ошибка: введено некорректное значение давления." << endl;
            return 0;
        }
        force.setForce(p);
        std::vector<ForceCondition3DPointer> boundaryForces;
        boundaryForces.push_back(&force);

        if (task == 1)
        {
            fem = new HexahedralFEM (&mesh, params, boundaryForces, boundaryConditions);
        }
        else
        {
            std::vector<MechanicalParameters3D> layers;
            layers.push_back(params);
            layers.push_back(paramsInner);
            layers.push_back(params);

            fem = new HexahedralFEM (&mesh, layers, boundaryForces, boundaryConditions); // многослойный расчет
        }

    }
    else if (task == 3)
    {   // spacecraft
        std::vector<FEMCondition3DPointer> boundaryConditions;
        SpaceBoundaryX sbx;
        SpaceBoundaryY sby;
        SpaceBoundaryZ sbz;
        boundaryConditions.push_back(&sbx);
        boundaryConditions.push_back(&sby);
        boundaryConditions.push_back(&sbz);
        std::vector<ForceCondition3DPointer> forces;
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

        fem = new HexahedralFEM (&mesh, layers, forces, boundaryConditions); // многослойный расчет
    }
    else if (task == 4)
    {
        MechanicalParameters3D params(65000.0, 0.3);
        TankBoundaryFixed fixed;
        TankBoundaryXMove xMove;
        std::vector<FEMCondition3DPointer> boundaryConditions;
        std::vector<ForceCondition3DPointer> boundaryForces;
        TankForce tankForce;

        int boundary_type = 0;
        cout << "Выберите тип граничных условий (1 - защемеление; 2 - допускается движение вдоль образующей): ";
        cin >> boundary_type;
        if (boundary_type != 1 && boundary_type != 2)
        {
            cout << "Ошибка: некорректный код типа граничных условий." << endl;
            return 0;
        }
        if (boundary_type == 1)
            boundaryConditions.push_back(&fixed);
        else
            boundaryConditions.push_back(&xMove);

        double p = 0;
        cout << "Введите значение внутреннего давления (p > 0, атмостферы): ";
        cin >> p;
        if (p <= 0)
        {
            cout << "Ошибка: введено некорректное значение давления." << endl;
            return 0;
        }
        tankForce.setPressure(p);
        boundaryForces.push_back(&tankForce);

        fem = new HexahedralFEM (&mesh, params, boundaryForces, boundaryConditions);
    }
    else if (task == 5)
    {
//        double strn[] =  {0.002,    0.0024, 0.003,  0.004,  0.0055, 0.0079, 0.013,  0.015,  0.0248, 0.032,  0.0361, 0.05,   0.1,    0.15};
//        std::vector<double> strain(strn, strn + sizeof(strn) / sizeof(strn[0]));
//        double strs[] =  {130,      140,    150,    160,    170,    180,    190,    200,    210,    220,    230,    270,    309,    330};
//        std::vector<double> stress(strs, strs + sizeof(strs) / sizeof(strs[0]));
        double strn[] =  {0.0026,    0.003, 0.0032, 0.0033,  0.0036,    0.0038, 0.0041, 0.0052, 0.0067, 0.0079, 0.01, 0.016,   0.026,   0.032,  0.11};
        std::vector<double> strain(strn, strn + sizeof(strn) / sizeof(strn[0]));
        double strs[] =  {180,      200,    210,    220,    240,        250,    260,    280,    300,    310,    320,    332,    350,    360,    400};
        std::vector<double> stress(strs, strs + sizeof(strs) / sizeof(strs[0]));
        TankBoundaryFixed fixed;
        TankBoundaryXMove xMove;
        std::vector<FEMCondition3DPointer> boundaryConditions;
        std::vector<ForceCondition3DPointer> boundaryForces;
        TankForce tankForce;

        int boundary_type = 0;
        cout << "Выберите тип граничных условий (1 - защемеление; 2 - допускается движение вдоль образующей): ";
        cin >> boundary_type;
        if (boundary_type != 1 && boundary_type != 2)
        {
            cout << "Ошибка: некорректный код типа граничных условий." << endl;
            return 0;
        }
        if (boundary_type == 1)
            boundaryConditions.push_back(&fixed);
        else
            boundaryConditions.push_back(&xMove);

        double p = 0;
        cout << "Введите значение внутреннего давления (p > 0, атмостферы): ";
        cin >> p;
        if (p <= 0)
        {
            cout << "Ошибка: введено некорректное значение давления." << endl;
            return 0;
        }
        tankForce.setPressure(p);
        boundaryForces.push_back(&tankForce);

        PlasticFem pFem(&mesh, strain, stress, 0.3, boundaryForces, boundaryConditions);

        cout << "Сохранение результатов в файл" << endl;
        save_vector(pFem.u(), "u.txt");
        save_vector(pFem.v(), "v.txt");
        save_vector(pFem.w(), "w.txt");
        save_vector(pFem.sigmaX(), "sigma_x.txt");
        save_vector(pFem.sigmaY(), "sigma_y.txt");
        save_vector(pFem.sigmaZ(), "sigma_z.txt");
        save_vector(pFem.tauXY(), "tau_xy.txt");
        save_vector(pFem.tauYZ(), "tau_yz.txt");
        save_vector(pFem.tauZX(), "tau_zx.txt");
        save_vector(pFem.sigma(), "sigma.txt");
        // сохране зон пластичности в файл
        fstream out;
        out.open("plastic_zones.txt", fstream::out);
        out << mesh.elementsCount() << endl;
        for (UInteger i = 0; i < mesh.elementsCount(); i++)
        {
            out << mesh.layer(i) << endl;
        }
        out.close();
    }

    if (fem)
    {
        cout << "Сохранение результатов в файл" << endl;
        save_vector(fem->u(), "u.txt");
        save_vector(fem->v(), "v.txt");
        save_vector(fem->w(), "w.txt");
        save_vector(fem->sigmaX(), "sigma_x.txt");
        save_vector(fem->sigmaY(), "sigma_y.txt");
        save_vector(fem->sigmaZ(), "sigma_z.txt");
        save_vector(fem->tauXY(), "tau_xy.txt");
        save_vector(fem->tauYZ(), "tau_yz.txt");
        save_vector(fem->tauZX(), "tau_zx.txt");
        save_vector(fem->sigma(), "sigma.txt");

        delete fem;
    }

    return 0;
}
