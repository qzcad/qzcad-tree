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
#include <functional>
#include <list>

#include "mesh.h"
using namespace msh;

#include "doublevector.h"
#include "doublematrix.h"
#include "mappeddoublematrix.h"
#include "femcondition.h"
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
    Fem(Mesh *mesh, UInteger freedom, const std::list<FemCondition *> &conditions);
    /**
     * @brief Деструктор
     */
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
     * @brief Метод строит конечно-элементное решение (линейное по умолчанию)
     */
    virtual void solve();
    /**
     * @brief Сигнатура функции для фильтрации узлов при формировании отчета
     */
    typedef std::function<bool(PointPointer)> PointFilterFunc;
    /**
     * @brief метод печатает в стандартный вывод результат отбора точек по значений в точках по указанному критерию
     * @param func Функция для отбора узлов
     */
    void reportNodeValues(PointFilterFunc func);
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
     * @brief Процедура ансамля локальной матрицы в глобальную
     * @param element Указатель на элемент сетки
     * @param local Локальная матрица
     * @param global Ссылка на глобальную матрицу
     */
    void assembly(ElementPointer element, const DoubleMatrix &local);
    /**
     * @brief Метод для учета граничных условий (начальных перемещений, температур и т.п.)
     * @param rowNumber Номер строки, которой соответствует начальное значение
     * @param value Начальное значение
     */
    void setInitialNodalValue(const UInteger &rowNumber, const double &value);
    /**
     * @brief Метод обработки списка граничных условий
     */
    void processInitialValues();
    /**
     * @brief Метод для решения СЛАУ в МКЭ
     */
    DoubleVector solveLinearSystem(bool cg=true);
    /**
     * @brief Метод для построения глобальной матрицы системы
     */
    virtual void buildGlobalMatrix() = 0;
    /**
     * @brief Метод для построения вектора системы
     */
    virtual void buildGlobalVector() = 0;
    /**
     * @brief Метод для обработки результтов решения
     * @param nodalValues
     */
    virtual void processSolution(const DoubleVector &nodalValues) = 0;

protected:
    Mesh *mesh_; //!< Указатель на сетку конечных элементов
    UInteger freedom_; //!< Количество степеней свободы
    UInteger dimension_; //!< Размерность задачи
    MappedDoubleMatrix global_; //!< Глобальная матрица жесткости в сжатом виде
    DoubleVector force_; //!< Вектор сил
    std::list<FemCondition *> conditions_; //!< Список граничныйх условий
};

#endif // FEM_H
