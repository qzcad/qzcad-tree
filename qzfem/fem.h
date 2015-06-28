/**
  * @author Сергей Чопоров
  * @date 17/06/2015
  * @version 1.0.0
  * @copyright Copyright 2015 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  * */
#ifndef FEM_H
#define FEM_H

#include <string>
#include <vector>

#include "mesh.h"
using namespace msh;

#include "doublevector.h"
#include "mappeddoublematrix.h"
using namespace mtx;

/**
 * @brief Базовый класс для реализаций МКЭ
 */
class Fem
{
public:
    /**
     * @brief Конструктор
     * @param mesh Указатель на сетку конечных элементов
     */
    Fem(Mesh *mesh);
    virtual ~Fem(){}
    /**
     * @brief Метод возвращает указатель на сетку конечных элементов
     * @return Указатель на сетку конечных элементов
     */
    Mesh *mesh();
    /**
     * @brief Метод возвращает количество степеней свободы
     * Количество степеней свободы определяет число векторов узловых значений (по одному вектору на степень свободы)
     * @return Количество степеней свободы
     */
    UInteger freedom() const;
    UInteger nodesVectorsCount() const;
    /**
     * @brief Метод возвращает узловые значения вектора с указанным номером
     * @param num Номер вектора узловых значений
     * @return Вектор значений
     */
    std::vector<double> nodeVector(UInteger num) const;
    /**
     * @brief Метод возвращает название вектора узловых значений
     * @param num Номер вектора узловых значений
     * @return Название вектора узловых значений
     */
    std::string nodeVectorName(UInteger num) const;
    /**
     * @brief Метод печатает в стандартный вывод экстремальные значения векторов узловых значений
     */
    void printNodeValuesExtremums() const;
protected:
    /**
     * @brief Метод для генерации квадратур для интегрирования на отрезке [-1; 1]
     * @param count Количество точек в квадратуре (от 1 до 5)
     * @param point Массив для записи координат на отрезке
     * @param weight Массив для записи весовых коэффициентов соответствующих координат
     */
    void quadrature(int count, DoubleVector &point, DoubleVector &weight);
    /**
     * @brief Метод для генерации квадратур для интегрирования на отрезке [0; 1] для интегрирования на треугольнике
     * @param count Количество точек в квадратуре (1, 3, 4)
     * @param xi Координата первого направления
     * @param eta Координата второго направления
     * @param weight Массив для записи весовых коэффициентов соответствующих координат
     */
    void quadrature(int count, DoubleVector &xi, DoubleVector &eta, DoubleVector &weight);
    /**
     * @brief Метод для учета граничных условий (начальных перемещений, температур и т.п.)
     * @param global Ссылка на глобальную матрицу жесткости
     * @param force Ссылка на вектор-столбец (сила, температура и т.п.)
     * @param rowNumber Номер строки, которой соответствует начальное значение
     * @param value Начальное значение
     */
    void setInitialNodalValue(MappedDoubleMatrix &global, DoubleVector &force, const UInteger &rowNumber, const double value);
    /**
     * @brief Метод для решения СЛАУ в МКЭ
     * @param global Ссылка на глобальную матрицу жесткости
     * @param force Ссылка на вектор-столбец
     */
    DoubleVector solve(MappedDoubleMatrix &global, DoubleVector &force);
protected:
    Mesh *mesh_; //!< Указатель на сетку конечных элементов
    UInteger freedom_; //!< Количество степеней свободы
    class NamedVector
    {
    public:
        NamedVector(const std::string &n, const std::vector<double> &v)
        {
            name = n;
            values = v;
        }

        std::string name;
        std::vector<double> values;
    };
    std::vector<NamedVector> nodeValues_;
};

#endif // FEM_H
