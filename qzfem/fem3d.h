#ifndef FEM3D_H
#define FEM3D_H

#include "fem.h"

using namespace msh;

#include "doublevector.h"
using namespace mtx;
/**
 * @brief Базовый класс для решения трехмерных задач МКЭ
 */
class Fem3D : public Fem
{
public:
    Fem3D(Mesh *mesh, UInteger freedom, const std::list<FemCondition *> &conditions);
protected:
    /**
     * @brief Метод для построения значений функций формы билинейного шестигранного элемента
     * @param xi Значение параметра первого направления местной системы координат
     * @param eta Значения параметра второго направления местной системы координат
     * @param x Массив x-координат узлов
     * @param y Массив y-координат узлов
     * @param N Значения функций формы
     * @param dNdX Значения x-производной функций формы
     * @param dNdY Значения y-производной функций формы
     * @return Якобиан преобразования в местную систему координат
     */
    static double isoHex8(const double &xi, const double &eta, const double &mu, const DoubleVector &x, const DoubleVector &y, const DoubleVector &z, DoubleVector &N, DoubleVector &dNdX, DoubleVector &dNdY, DoubleVector &dNdZ);
    /**
     * @brief Метод для построения значений функций формы билинейного тетраэдричесого элемента
     * @param xi Значение параметра первого направления местной системы координат
     * @param eta Значения параметра второго направления местной системы координат
     * @param x Массив x-координат узлов
     * @param y Массив y-координат узлов
     * @param N Значения функций формы
     * @param dNdX Значения x-производной функций формы
     * @param dNdY Значения y-производной функций формы
     * @return Якобиан преобразования в местную систему координат
     */
    static double isoTet4(const double &xi, const double &eta, const double &mu, const DoubleVector &x, const DoubleVector &y, const DoubleVector &z, DoubleVector &N, DoubleVector &dNdX, DoubleVector &dNdY, DoubleVector &dNdZ);
};

#endif // FEM3D_H
