#include <iostream>
#include "plasticfem.h"

PlasticFem::PlasticFem(HexahedralMesh3D* mesh, const std::vector<double> &strain, const std::vector<double> &stress, const double &nu, const std::vector<FEMCondition3DPointer> &boundaryForces, const std::vector<FEMCondition3DPointer> &boundaryConditions)
{
    std::vector<MechanicalParameters3D> layers;

    for (int i = 0; i < stress.size(); i++)
    {
        E[i] = stress[i] / strain[i];
        layers.push_back(MechanicalParameters3D(stress[i] / strain[i], nu));
    }
    HexahedralFEM *hFem;
    std::cout << "Нулевое приближение..." << std::endl;
    hFem = new HexahedralFEM(mesh, layers[0], boundaryForces,boundaryConditions);
    std::vector<double> sigma = hFem->sigma();
    mesh->clearElementValues();
    for (msh::UInteger i = 0; i < mesh->elementsCount(); i++) mesh->pushElementValue(0);
    for ()
}
