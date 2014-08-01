/**
  * Проект для консольного тестирования конечно элементных расчетов
  * */
#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/numeric/ublas/io.hpp>
#include "hexahedralmesh3d.h"
#include "floatingmatrix.h"

using namespace std;
using namespace msh;

int main()
{
    cout << "Программа для тестирования расчетов при помощи МКЭ" << endl;
    fstream input;
    UInteger nodesCount;
    int freedom = 0;
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
        UInteger p[8];
        for (int j = 0; j < 8; j++)
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
    cout << "Внешние слоя: E = " << outerE << "; nu = " << outerNu << "." << endl;



    return 0;
}

