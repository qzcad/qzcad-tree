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
    virtual std::string nodeVectorName(UInteger num) const = 0;
    /**
     * @brief Метод возвращает количество векторов значений, определенных на элементе
     * @return Количество векторов значений, определенных на элементе
     */
    UInteger elementVectorsCount() const;
    /**
     * @brief Метод возвращает значения, определенные на элементе, для вектора с указанным номером
     * @param num Номер вектора значений, определенных на элементе
     * @return Вектор значений
     */
    std::vector<double> elementVector(UInteger num) const;
    /**
     * @brief Метод возвращает название вектора значений, определенных на элементе
     * @param num Номер вектора значений, определеных на элементе
     * @return Название вектора значений, определенных на элеменете
     */
    virtual std::string elementVectorName(UInteger num) const = 0;
    /**
     * @brief Метод печатает в стандартный вывод экстремальные значения векторов узловых значений
     */
    void printNodeValuesExtremums() const;
    /**
     * @brief Метод печатает в стандартный вывод экстремальные значения векторов значений, определенных на элементе
     */
    void printElementValuesExtremums() const;
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
    void solve(MappedDoubleMatrix &global, DoubleVector &force);
protected:
    Mesh *mesh_; //!< Указатель на сетку конечных элементов
    UInteger freedom_; //!< Количество степеней свободы (= количество векторов узловых значений)
    DoubleVector nodeValues_; //!< Вектор узловых значений (размер freedom * nodesCount)
    UInteger elementVectorsCount_; //!< Количество векторов значений, определенных на элементе
    DoubleVector elementValues_; //!< Вектор значений, определенных на элементе (размер elementVectorsCount_ * elementsCount)
};

#endif // FEM_H
