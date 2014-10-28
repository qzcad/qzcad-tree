#include <iostream>
#include "plasticfem.h"

PlasticFem::PlasticFem(HexahedralMesh3D* mesh, const std::vector<double> &strain, const std::vector<double> &stress, const double &nu, const std::vector<ForceCondition3DPointer> &forceCondition, const std::vector<FEMCondition3DPointer> &boundaryConditions)
{
    const int stress_size = stress.size();
    std::vector<MechanicalParameters3D> layers;
    int it = 0; // начение множителя
    bool isRupture = false;

    for (int i = 0; i < stress_size; i++)
    {
        layers.push_back(MechanicalParameters3D(stress[i] / strain[i], nu));
    }
    HexahedralFEM *hFem;
    std::cout << "Нулевое приближение..." << std::endl;
    hFem = new HexahedralFEM(mesh, layers[0], forceCondition, boundaryConditions);
    sigma_ = hFem->sigma();
    sigmaX_ = hFem->sigmaX();
    sigmaY_ = hFem->sigmaY();
    sigmaZ_ = hFem->sigmaZ();
    tauXY_ = hFem->tauXY();
    tauYZ_ = hFem->tauYZ();
    tauZX_ = hFem->tauZX();
    u_ = hFem->u();
    v_ = hFem->v();
    w_ = hFem->w();
    delete hFem;
    std::cout << "Нулевое приближение построено." << std::endl;
    // @ цель этого участка кода "перепрыгнуть" зону упругости (линейный участок)
    double sigmaMax = 0.0;
    for (msh::UInteger i = 0; i < mesh->elementsCount(); i++)
    {
        if (sigmaMax < sigma_[i]) sigmaMax = sigma_[i];
    }
    it = (int)(stress[0] / sigmaMax); // во сколько раз начальная интенсивность меньше максимума для зоны упругости
    std::cout << "Множитель линейного участка: " << it << std::endl;
    vectorScale(sigma_, (double)it);
    vectorScale(sigmaX_, (double)it);
    vectorScale(sigmaY_, (double)it);
    vectorScale(sigmaZ_, (double)it);
    vectorScale(tauXY_, (double)it);
    vectorScale(tauYZ_, (double)it);
    vectorScale(tauZX_, (double)it);
    vectorScale(u_, (double)it);
    vectorScale(v_, (double)it);
    vectorScale(w_, (double)it);
    // @
    mesh->clearLayers();
    for (msh::UInteger i = 0; i < mesh->elementsCount(); i++) mesh->pushLayer(0);
    while (!isRupture)
    {
        for (typename std::vector<double>::size_type i = 0; i < sigma_.size(); i++)
        {
            if (stress[stress_size - 1] <= sigma_[i])
            {
                isRupture = true;
                mesh->setLayer(i, mesh->layer(i) + 1);
            }
            else
            {
                for (int j = 0; j < stress_size - 1; j++)
                {
                    if (stress[j] < sigma_[i])
                    {
                        mesh->setLayer(i, j + 1);
                    }
                }
            }
        }
        if (!isRupture)
        {
            ++it;
            hFem = new HexahedralFEM(mesh, layers, forceCondition, boundaryConditions);
            vectorSumAB(u_, hFem->u());
            vectorSumAB(v_, hFem->v());
            vectorSumAB(w_, hFem->w());
            vectorSumAB(sigmaX_, hFem->sigmaX());
            vectorSumAB(sigmaY_, hFem->sigmaY());
            vectorSumAB(sigmaZ_, hFem->sigmaZ());
            vectorSumAB(tauXY_, hFem->tauXY());
            vectorSumAB(tauYZ_, hFem->tauYZ());
            vectorSumAB(tauZX_, hFem->tauZX());
            vectorSumAB(sigma_, hFem->sigma());
            delete hFem;

            std::cout << "Итерация " << it + 1 << " завершена. Множитель: " << it << std::endl;
        }
    }
    std::cout << "Конечно-элементный анализ завершен. Множитель: " << it << std::endl;
}

void PlasticFem::vectorSumAB(std::vector<double> &a, const std::vector<double> &b)
{
    for (typename std::vector<double>::size_type i = 0; i < a.size(); i++)
    {
        a[i] += b[i];
    }
}

void PlasticFem::vectorScale(std::vector<double> &a, double b)
{
    for (typename std::vector<double>::size_type i = 0; i < a.size(); i++)
    {
        a[i] *= b;
    }
}

std::vector<double> PlasticFem::sigma() const
{
    return sigma_;
}

void PlasticFem::setSigma(const std::vector<double> &sigma)
{
    sigma_ = sigma;
}

std::vector<double> PlasticFem::tauZX() const
{
    return tauZX_;
}

void PlasticFem::setTauZX(const std::vector<double> &tauZX)
{
    tauZX_ = tauZX;
}

std::vector<double> PlasticFem::tauYZ() const
{
    return tauYZ_;
}

void PlasticFem::setTauYZ(const std::vector<double> &tauYZ)
{
    tauYZ_ = tauYZ;
}

std::vector<double> PlasticFem::tauXY() const
{
    return tauXY_;
}

void PlasticFem::setTauXY(const std::vector<double> &tauXY)
{
    tauXY_ = tauXY;
}

std::vector<double> PlasticFem::sigmaZ() const
{
    return sigmaZ_;
}

void PlasticFem::setSigmaZ(const std::vector<double> &sigmaZ)
{
    sigmaZ_ = sigmaZ;
}

std::vector<double> PlasticFem::sigmaY() const
{
    return sigmaY_;
}

void PlasticFem::setSigmaY(const std::vector<double> &sigmaY)
{
    sigmaY_ = sigmaY;
}

std::vector<double> PlasticFem::sigmaX() const
{
    return sigmaX_;
}

void PlasticFem::setSigmaX(const std::vector<double> &sigmaX)
{
    sigmaX_ = sigmaX;
}

std::vector<double> PlasticFem::w() const
{
    return w_;
}

void PlasticFem::setW(const std::vector<double> &w)
{
    w_ = w;
}

std::vector<double> PlasticFem::v() const
{
    return v_;
}

void PlasticFem::setV(const std::vector<double> &v)
{
    v_ = v;
}

std::vector<double> PlasticFem::u() const
{
    return u_;
}

void PlasticFem::setU(const std::vector<double> &u)
{
    u_ = u;
}

