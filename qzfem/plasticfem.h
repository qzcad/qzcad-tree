/**
  * @author Сергей Чопоров
  * @date 06/10/2014
  * @version 1.0.0
  * */
#ifndef PLASTICFEM_H
#define PLASTICFEM_H
#include <vector>
#include "hexahedralfem.h"
/**
 * @brief Класс для расчета конечно-элементных моделей в упруго-пластичной/пластичной постановке
 */
class PlasticFem
{
public:
    /**
     * @brief Конструктор для упруго-пластичного расчета. Шастигранные конечные элементы
     * @param mesh Указатель на сетку шестигранных элементов
     * @param strain Массив деформаций
     * @param stress Соответствующий массиву деформаций массив напряжений
     * @param nu Коэффициент Пуассона
     * @param forceCondition Массив сил (шаг)
     * @param boundaryConditions Массив граничных условий
     * Номера слоев будут соответствовать зонам пластичности
     */
    PlasticFem(HexahedralMesh3D* mesh, const std::vector<double> &strain, const std::vector<double> &stress, const double &nu, const std::vector<ForceCondition3DPointer> &forceCondition, const std::vector<FEMCondition3DPointer> &boundaryConditions);
    /// @brief getters and setters
    /// @{
    std::vector<double> u() const;
    void setU(const std::vector<double> &u);

    std::vector<double> v() const;
    void setV(const std::vector<double> &v);

    std::vector<double> w() const;
    void setW(const std::vector<double> &w);

    std::vector<double> sigmaX() const;
    void setSigmaX(const std::vector<double> &sigmaX);

    std::vector<double> sigmaY() const;
    void setSigmaY(const std::vector<double> &sigmaY);

    std::vector<double> sigmaZ() const;
    void setSigmaZ(const std::vector<double> &sigmaZ);

    std::vector<double> tauXY() const;
    void setTauXY(const std::vector<double> &tauXY);

    std::vector<double> tauYZ() const;
    void setTauYZ(const std::vector<double> &tauYZ);

    std::vector<double> tauZX() const;
    void setTauZX(const std::vector<double> &tauZX);

    std::vector<double> sigma() const;
    void setSigma(const std::vector<double> &sigma);
    /// @}
protected:
    /**
     * @brief a = a + b
     * @param a Выходной параметр (передача по ссылке)
     * @param b Входной параметр
     */
    void vectorSumAB(std::vector<double> &a, const std::vector<double> &b);
    /**
     * @brief Увеличить значение каждого элемента вектора a в b раз
     * @param a Вектор (входной-выходной параметр)
     * @param b Множитель
     */
    void vectorScale(std::vector<double> &a, double b);
private:
    std::vector<double> u_; //!< Перемещения в первом направлении (x)
    std::vector<double> v_; //!< Перемещения во втором направлении (y)
    std::vector<double> w_; //!< Перемещения в третьем направлении (z)
    std::vector<double> sigmaX_; //!< Нормальные компоненты напряжения, параллельные первому направлению (x)
    std::vector<double> sigmaY_; //!< Нормальные компоненты напряжения, параллельные второму направлению (y)
    std::vector<double> sigmaZ_; //!< Нормальные компоненты напряжения, параллельные третьему направлению (z)
    std::vector<double> tauXY_; //!< Касательные компоненты напряжения, в плоскости 1-2 (x-y)
    std::vector<double> tauYZ_; //!< Касательные компоненты напряжения, в плоскости 2-3 (y-z)
    std::vector<double> tauZX_; //!< Касательные компоненты напряжения, в плоскости 3-1 (z-x)
    std::vector<double> sigma_; //!< Интенсивность напряжений
};

#endif // PLASTICFEM_H
