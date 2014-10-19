#include <iostream>
#include "plasticfem.h"

PlasticFem::PlasticFem(HexahedralMesh3D* mesh, const std::vector<double> &strain, const std::vector<double> &stress, const double &nu, const std::vector<FEMCondition3DPointer> &boundaryForces, const std::vector<FEMCondition3DPointer> &boundaryConditions)
{
    const int stress_size = stress.size();
    std::vector<MechanicalParameters3D> layers;
    int it = 1;
    bool isRupture = false;

    for (int i = 0; i < stress_size; i++)
    {
        layers.push_back(MechanicalParameters3D(stress[i] / strain[i], nu));
    }
    HexahedralFEM *hFem;
    std::cout << "Нулевое приближение..." << std::endl;
    hFem = new HexahedralFEM(mesh, layers[0], boundaryForces, boundaryConditions);
    std::vector<double> sigma = hFem->sigma();
    delete hFem;
    mesh->clearLayers();
    for (msh::UInteger i = 0; i < mesh->elementsCount(); i++) mesh->pushLayer(0);
    while (!isRupture)
    {
        for (typename std::vector<double>::size_type i = 0; i < sigma.size(); i++)
        {
            if (stress[stress_size - 1] <= sigma[i])
            {
                isRupture = true;
                mesh->setLayer(i, mesh->layer(i) + 1);
            }
            else
            {
                for (int j = 0; j < stress_size - 1; j++)
                {
                    if (stress[j] < sigma[i])
                    {
                        mesh->setLayer(i, j + 1);
                    }
                }
            }
        }
        hFem = new HexahedralFEM(mesh, layers, boundaryForces, boundaryConditions);
        std::vector<double> sigma_new = hFem->sigma();
        delete hFem;
        for (typename std::vector<double>::size_type i = 0; i < sigma.size(); i++)
        {
            sigma[i] += sigma_new[i];
        }
        std::cout << "Итерация " << it << " завершена." << std::endl;
        it++;
    }
}
