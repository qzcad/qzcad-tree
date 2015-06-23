#ifndef FEM2D_H
#define FEM2D_H

#include <functional>

#include "fem.h"

#include "point2d.h"
using namespace msh;

#include "doublevector.h"
using namespace mtx;
/**
 * @brief Базовый класс для решения двумерных задач МКЭ
 */
class Fem2D : public Fem
{
public:
    /**
     * @brief Конструктор
     * @param mesh Указатель на сетку конечных элементов
     */
    Fem2D(Mesh *mesh);
    /**
     * @brief Сигнатура функции, описывающей направление действия начального условия: -1 - словие не применять; 0 - применять во всех направления свободы; 1 - только первое направление и т.д.
     */
    typedef std::function<int(double, double)> BoundaryConditionFunction;
    /**
     * @brief Сигнатура функции, описывающей действие вектора в точке
     */
    typedef std::function<Point2D(double, double)> VectorFunction2D;
protected:
    /**
     * @brief Метод для построения значений функций формы билинейного четырехугольного элемента
     * @param xi Значение параметра первого направления местной системы координат
     * @param eta Значения параметра второго направления местной системы координат
     * @param x Массив x-координат узлов
     * @param y Массив y-координат узлов
     * @param N Значения функций формы
     * @param dNdX Значения x-производной функций формы
     * @param dNdY Значения y-производной функций формы
     * @return Якобиан преобразования в местную систему координат
     */
    double isoQuad4(const double &xi, const double &eta, double x[], double y[],
                    DoubleVector &N, DoubleVector &dNdX, DoubleVector &dNdY);
    /**
     * @brief Метод для построения значений функций формы линейного треугольного элемента
     * @param xi Значение параметра первого направления местной системы координат
     * @param eta Значения параметра второго направления местной системы координат
     * @param x Массив x-координат узлов
     * @param y Массив y-координат узлов
     * @param N Значения функций формы
     * @param dNdX Значения x-производной функций формы
     * @param dNdY Значения y-производной функций формы
     * @return Якобиан преобразования в местную систему координат
     */
    double isoTriangle3(const double &xi, const double &eta, double x[], double y[],
                        DoubleVector &N, DoubleVector &dNdX, DoubleVector &dNdY);
};

#endif // FEM2D_H
