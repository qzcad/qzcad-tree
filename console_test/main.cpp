/**
  * Проект для консольного тестирования конечно элементных расчетов
  * */
#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/numeric/ublas/io.hpp>
#include <boost/progress.hpp>
#include "hexahedralmesh3d.h"
#include "floatingmatrix.h"
#include "globalmatrix.h"
#include "floatingvector.h"

using namespace std;
using namespace msh;

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
    }

    return 0;
}

