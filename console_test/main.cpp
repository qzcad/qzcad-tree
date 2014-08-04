/**
  * Проект для консольного тестирования конечно элементных расчетов
  * */
#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/progress.hpp>
#include "hexahedralmesh3d.h"
#include "floatingmatrix.h"
#include "globalmatrix.h"
#include "floatingvector.h"
#include "invertmatrix.hpp"
#include "determinant.hpp"
#include "conjugategradient.hpp"

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

int main()
{
    cout << "Программа для тестирования расчетов при помощи МКЭ" << endl;
    fstream input;
    UInteger nodesCount;
    int freedom = 0;
    const int elementNodes = 8;
    HexahedralMesh3D mesh;
    UInteger elementsCount;

    const Floating innerG = 2.77E+4;
    const Floating innerK = 6.0E+4;
    const Floating innerNu = (3.0 * innerK - 2.0 * innerG) / (2.0 * (3.0 * innerK + innerG)); // коэффициент Пуассона среднего слоя
    const Floating innerE = 2.0 * innerG * (1.0 + innerNu); // модуль Юнга среднего слоя
    FloatingMatrix innerD(6, 6);

    const Floating outerG = 8.0E+4;
    const Floating outerNu = 0.27; // коэффициент Пуассона внешних слоев
    const Floating outerE = 2.0 * outerG * (1.0 + outerNu); // модуль Юнга внешних слоев
    FloatingMatrix outerD(6, 6);

    const Floating H = 0.016; // толщина заполнителя (внутреннего слоя)
    const Floating delta = 0.001; // толщина внешних слоев
    const Floating outerR = 0.4;
    const Floating innerR = 0.2;

    const int gaussCount = 2; // количество точек в квадратурах для интегрирования
    Floating gaussPoint[] = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)}; // координаты точек квадратуры
    Floating gaussWeight[] = {1.0, 1.0}; // веса точек квадратуры

    const Floating q = -0.09; // интенсивность нагрузки

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

    innerD(0, 0) = 1.0; innerD(0, 1) = innerNu / (1.0 - innerNu); innerD(0, 2) = innerNu / (1.0 - innerNu);
    innerD(1, 0) = innerNu / (1.0 - innerNu); innerD(1, 1) = 1.0; innerD(1, 2) = innerNu / (1.0 - innerNu);
    innerD(2, 0) = innerNu / (1.0 - innerNu); innerD(2, 1) = innerNu / (1.0 - innerNu); innerD(2, 2) = 1.0;
    innerD(3, 3) = (1.0 - 2.0 * innerNu) / (2.0 * (1.0 - innerNu));
    innerD(4, 4) = (1.0 - 2.0 * innerNu) / (2.0 * (1.0 - innerNu));
    innerD(5, 5) = (1.0 - 2.0 * innerNu) / (2.0 * (1.0 - innerNu));
    innerD = innerE * (1.0 - innerNu) / ((1.0 + innerNu) * (1.0 - 2.0 * innerNu)) * innerD;

    cout << "Средний слой: E = " << innerE << "; nu = " << innerNu << ";" << endl;
    cout << "D = " << endl << innerD << endl;

    outerD(0, 0) = 1.0; outerD(0, 1) = outerNu / (1.0 - outerNu); outerD(0, 2) = outerNu / (1.0 - outerNu);
    outerD(1, 0) = outerNu / (1.0 - outerNu); outerD(1, 1) = 1.0; outerD(1, 2) = outerNu / (1.0 - outerNu);
    outerD(2, 0) = outerNu / (1.0 - outerNu); outerD(2, 1) = outerNu / (1.0 - outerNu); outerD(2, 2) = 1.0;
    outerD(3, 3) = (1.0 - 2.0 * outerNu) / (2.0 * (1.0 - outerNu));
    outerD(4, 4) = (1.0 - 2.0 * outerNu) / (2.0 * (1.0 - outerNu));
    outerD(5, 5) = (1.0 - 2.0 * outerNu) / (2.0 * (1.0 - outerNu));
    outerD = outerE * (1.0 - outerNu) / ((1.0 + outerNu) * (1.0 - 2.0 * outerNu)) * outerD;

    cout << "Внешние слоя: E = " << outerE << "; nu = " << outerNu << ";" << endl;
    cout << "D = " << endl << outerD << endl;

    cout << "Построение глобальной матрицы..." << endl;

    UInteger systemDimension = nodesCount * freedom;
    GlobalMatrix globalMatrix (systemDimension, systemDimension);
    FloatingVector force(systemDimension);
    FloatingVector displacement(systemDimension);

    boost::progress_display progressBar(elementsCount);
    for (UInteger elementNumber = 0; elementNumber < elementsCount; elementNumber++)
    {
        ++progressBar;

        Floating x[elementNodes];
        Floating y[elementNodes];
        Floating z[elementNodes];
        FloatingMatrix localMatrix(elementNodes * freedom, elementNodes * freedom);
        UInteger index_i = 0;
        UInteger index_j = 0;
        for (int i = 0; i < elementNodes * freedom; i++)
        {
            for (int j = 0; j < elementNodes * freedom; j++)
            {
                localMatrix(i, j) = 0.0;
            }
        }
        // извлечение координат узлов
        ElementPointer element = mesh.element(elementNumber);
        for (int i = 0; i < elementNodes; i++)
        {
            PointPointer point = mesh.node(element->vertexNode(i));
            x[i] = point->x();
            y[i] = point->y();
            z[i] = point->z();
        }
        // интегрирование квадратурами
        for (int iXi = 0; iXi < gaussCount; iXi++)
        {
            Floating xi = gaussPoint[iXi];
            Floating xiWeight = gaussWeight[iXi];
            for (int iEta = 0; iEta < gaussCount; iEta++)
            {
                Floating eta = gaussPoint[iEta];
                Floating etaWeight = gaussWeight[iEta];
                for (int iMu = 0; iMu < gaussCount; iMu++)
                {
                    Floating mu = gaussPoint[iMu];
                    Floating muWeight = gaussWeight[iMu];
                    // значения функций формы
                    FloatingVector N(elementNodes);
                    // значения производных функций формы в местных координатах
                    FloatingVector dNdXi(elementNodes);
                    FloatingVector dNdEta(elementNodes);
                    FloatingVector dNdMu(elementNodes);
                    // Якобиан
                    FloatingMatrix Jacobian(freedom, freedom);
                    FloatingMatrix invJacobian(freedom, freedom);
                    Floating detJacobian;
                    // значения производных функций формы в глобальных координатах
                    FloatingVector dNdX(elementNodes);
                    FloatingVector dNdY(elementNodes);
                    FloatingVector dNdZ(elementNodes);
                    // матрицы вариационной постановки
                    FloatingMatrix B(6, 24);
                    FloatingMatrix BT(24, 6);
                    FloatingMatrix BTD(24, 6);
                    short out = 0;

                    N(0) = (1 - xi) * (1 - eta) * (1 - mu) / 8.0;
                    N(1) = (1 - xi) * (1 - eta) * (1 + mu) / 8.0;
                    N(2) = (1 + xi) * (1 - eta) * (1 + mu) / 8.0;
                    N(3) = (1 + xi) * (1 - eta) * (1 - mu) / 8.0;
                    N(4) = (1 - xi) * (1 + eta) * (1 - mu) / 8.0;
                    N(5) = (1 - xi) * (1 + eta) * (1 + mu) / 8.0;
                    N(6) = (1 + xi) * (1 + eta) * (1 + mu) / 8.0;
                    N(7) = (1 + xi) * (1 + eta) * (1 - mu) / 8.0;

                    dNdXi(0) = -(1 - eta) * (1 - mu) / 8.0;
                    dNdXi(1) = -(1 - eta) * (1 + mu) / 8.0;
                    dNdXi(2) = (1 - eta) * (1 + mu) / 8.0;
                    dNdXi(3) = (1 - eta) * (1 - mu) / 8.0;
                    dNdXi(4) = -(1 + eta) * (1 - mu) / 8.0;
                    dNdXi(5) = -(1 + eta) * (1 + mu) / 8.0;
                    dNdXi(6) = (1 + eta) * (1 + mu) / 8.0;
                    dNdXi(7) = (1 + eta) * (1 - mu) / 8.0;

                    dNdEta(0) = -(1 - xi) * (1 - mu) / 8.0;
                    dNdEta(1) = -(1 - xi) * (1 + mu) / 8.0;
                    dNdEta(2) = -(1 + xi) * (1 + mu) / 8.0;
                    dNdEta(3) = -(1 + xi) * (1 - mu) / 8.0;
                    dNdEta(4) = (1 - xi) * (1 - mu) / 8.0;
                    dNdEta(5) = (1 - xi) * (1 + mu) / 8.0;
                    dNdEta(6) = (1 + xi) * (1 + mu) / 8.0;
                    dNdEta(7) = (1 + xi) * (1 - mu) / 8.0;

                    dNdMu(0) = -(1 - xi) * (1 - eta) / 8.0;
                    dNdMu(1) = (1 - xi) * (1 - eta) / 8.0;
                    dNdMu(2) = (1 + xi) * (1 - eta) / 8.0;
                    dNdMu(3) = -(1 + xi) * (1 - eta) / 8.0;
                    dNdMu(4) = -(1 - xi) * (1 + eta) / 8.0;
                    dNdMu(5) = (1 - xi) * (1 + eta) / 8.0;
                    dNdMu(6) = (1 + xi) * (1 + eta) / 8.0;
                    dNdMu(7) = -(1 + xi) * (1 + eta) / 8.0;

                    Jacobian(0, 0) = 0.0; Jacobian(0, 1) = 0.0; Jacobian(0, 2) = 0.0;
                    Jacobian(1, 0) = 0.0; Jacobian(1, 1) = 0.0; Jacobian(1, 2) = 0.0;
                    Jacobian(2, 0) = 0.0; Jacobian(2, 1) = 0.0; Jacobian(2, 2) = 0.0;

                    for (int i = 0; i < elementNodes; i++)
                    {
                        Jacobian(0, 0) += dNdXi(i) * x[i]; Jacobian(0, 1) += dNdXi(i) * y[i]; Jacobian(0, 2) += dNdXi(i) * z[i];
                        Jacobian(1, 0) += dNdEta(i) * x[i]; Jacobian(1, 1) += dNdEta(i) * y[i]; Jacobian(1, 2) += dNdEta(i) * z[i];
                        Jacobian(2, 0) += dNdMu(i) * x[i]; Jacobian(2, 1) += dNdMu(i) * y[i]; Jacobian(2, 2) += dNdMu(i) * z[i];
                    }
                    invertMatrix(Jacobian, invJacobian);
                    detJacobian = determinant(Jacobian);

                    for (int i = 0; i < elementNodes; i++)
                    {
                        dNdX(i) = invJacobian(0, 0) * dNdXi(i) + invJacobian(0, 1) * dNdEta(i) + invJacobian(0, 2) * dNdMu(i);
                        dNdY(i) = invJacobian(1, 0) * dNdXi(i) + invJacobian(1, 1) * dNdEta(i) + invJacobian(1, 2) * dNdMu(i);
                        dNdZ(i) = invJacobian(2, 0) * dNdXi(i) + invJacobian(2, 1) * dNdEta(i) + invJacobian(2, 2) * dNdMu(i);
                    }

                    for (int i = 0; i < 6; i++)
                    {
                        for (int j = 0; j < 24; j++)
                        {
                            B(i, j) = 0.0;
                        }
                    }
                    for (int i = 0; i < 8; i++)
                    {
                        B(0, i) = dNdX(i);
                        B(1, i + 8) = dNdY(i);
                        B(2, i + 16) = dNdZ(i);
                        B(3, i) = dNdY(i); B(3, i + 8) = dNdX(i);
                        B(4, i + 8) = dNdZ(i); B(4, i + 16) = dNdY(i);
                        B(5, i) = dNdZ(i); B(5, i + 16) = dNdX(i);
                    }
                    BT = trans(B);

                    for (int i = 0; i < elementNodes; i++)
                    {
                        if (y[i] < delta + 0.00001 || y[i] > H + delta - 0.00001) out++; // элемент находится во внешних слоях
                    }
                    if (out == 8)
                    {
                        BTD = prod(BT, outerD);
                    }
                    else
                    {
                        BTD = prod(BT, innerD);
                    }

                    localMatrix += prod(BTD, B) * detJacobian * xiWeight * etaWeight * muWeight;

                } // for iMu
            } // for iEta
        } // for iXi
        // Ансамблирование
        for (int i = 0; i < elementNodes * freedom; i++)
        {
            if (i < elementNodes)
            {
                index_i = element->vertexNode(i);
            }
            else if (i < 2 * elementNodes)
            {
                index_i = element->vertexNode(i - elementNodes) + nodesCount;
            }
            else
            {
                index_i = element->vertexNode(i - 2L * elementNodes) + 2L * nodesCount;
            }

            for (int j = 0; j < elementNodes * freedom; j++)
            {
                if (j < elementNodes)
                {
                    index_j = element->vertexNode(j);
                }
                else if (j < 2 * elementNodes)
                {
                    index_j = element->vertexNode(j - elementNodes) + nodesCount;
                }
                else
                {
                    index_j = element->vertexNode(j - 2L * elementNodes) + 2L * nodesCount;
                }
                globalMatrix(index_i, index_j) += localMatrix(i, j);
            } // for j
        } // for i
    } // for elementNumber
    cout << "Построение вектора нагрузок... " << endl;
    for (UInteger i = 0; i < systemDimension; i++)
        force(i) = 0.0;
    progressBar.restart(elementsCount);
    Floating area = 0.0;
    for (UInteger i = 0; i < elementsCount; i++)
    {
        ++progressBar;
        if (mesh.isBorderElement(i))
        {
            ElementPointer element = mesh.element(i);
            for (int j = 0; j < element->facesCount(); j++)
            {
                UIntegerVector face = element->face(j);
                bool isBorderFace = true;
                for (int k = 0; k < 4; k++)
                {
                    if (mesh.nodeType(face[k]) == INNER || !isEquil(mesh.node(face[k])->y(), H + 2.0 * delta))
                    {
                        isBorderFace = false;
                    }
                } // for k
                if (isBorderFace)
                {
                    Floating faceArea = mesh.faceArea(face);
                    area += faceArea;
                    for (int k = 0; k < 4; k++)
                    {
                        force(face[k] + nodesCount) = force(face[k] + nodesCount) + q * faceArea / 4.0;
                    } // for k
                }
            } // for j
        } // if
    } // for i
    cout << "Площадь нагруженной поверхноти: " << area << endl;
    cout << "Учет граничных условий..." << endl;
    progressBar.restart(nodesCount);
    for (UInteger i = 0; i < nodesCount; i++)
    {
        ++progressBar;
        PointPointer point = mesh.node(i);
#if FIXED == 1
        if ( isEquil(node[i].x * node[i].x + node[i].y * node[i].y, B * B) ) // защемление
        {
            for (long j = 0; j < system_dimension; j++)
            {
                if (GlobalMatrix(i, j) != 0.0) GlobalMatrix(i, j) = 0.0; // защемление
                if (GlobalMatrix(i + nodes_count, j) != 0.0) GlobalMatrix(i + nodes_count, j) = 0.0; // защемление
                if (GlobalMatrix(i + 2 * nodes_count, j) != 0.0) GlobalMatrix(i + 2 * nodes_count, j) = 0.0;
            }
            GlobalMatrix(i, i) = 1.0; // защемление
            force(i) = 0.0; // защемление
            GlobalMatrix(i + nodes_count, i + nodes_count) = 1.0; // защемление
            force(i + nodes_count) = 0.0; // защемление
            GlobalMatrix(i + 2 * nodes_count, i + 2 * nodes_count) = 1.0;
            force(i + 2 * nodes_count) = 0.0;
        }
#else
        if ( isEquil(point->y(), 0.0) && isEquil(point->x() * point->x() + point->z() * point->z(), outerR * outerR) ) // опирание
        {
            for (UInteger j = 0; j < systemDimension; j++)
            {
                if (globalMatrix(i + nodesCount, j) != 0.0) globalMatrix(i + nodesCount, j) = 0.0;
            } // j
            globalMatrix(i + nodesCount, i + nodesCount) = 1.0;
            force(i + nodesCount) = 0.0;
        }
#endif
        // для сектора цилиндра необходимо зафиксировать грани
        if ( isEquil(point->x(), 0.0) )
        {
            for (long j = 0; j < systemDimension; j++)
            {
                if (globalMatrix(i, j) != 0.0) globalMatrix(i, j) = 0.0;
            }
            globalMatrix(i, i) = 1.0;
            force(i) = 0.0;
        }
        if ( isEquil(point->z(), 0.0) )
        {
            for (long j = 0; j < systemDimension; j++)
            {
                if (globalMatrix(i + 2L * nodesCount, j) != 0.0) globalMatrix(i + 2L * nodesCount, j) = 0.0;
            }
            globalMatrix(i + 2L * nodesCount, i + 2L * nodesCount) = 1.0;
            force(i + 2L * nodesCount) = 0.0;
        }
    } // for i
    cout << "Ассоциация и частичная очистка памяти..." << endl;
    compressed_matrix<Floating> CGM(systemDimension, systemDimension);
    CGM.assign(globalMatrix);
    globalMatrix.clear();
    cout << "Решение СЛАУ..." << endl;
    int res = conjugateGradient(CGM, displacement,force,row_major_tag());
    cout << "Выполнено итераций в методу сопряженных градиентов: " << res;
    cout << "Обработка результатов..." << endl;
    progressBar.restart(systemDimension);
    Floating maxW = displacement(0);
    for (UInteger i = 0; i < systemDimension; i++)
    {
        ++progressBar;
        if (maxW < fabs(displacement(i)))
        {
            maxW = fabs(displacement(i));
        }
    } // for i
    cout << "Максимальное перемещение: " << maxW << endl;
    return 0;
}

