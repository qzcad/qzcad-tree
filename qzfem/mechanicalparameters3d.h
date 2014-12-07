/**
  * @author Сергей Чопоров
  * @date 14/08/2014
  * @version 1.0.1
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  * */
#ifndef MECHANICALPARAMETERS3D_H
#define MECHANICALPARAMETERS3D_H
/**
 * @brief Класс для передачи основных механических параметров (констант): модуль Юнга, коэффициент Пуассона
 */
class MechanicalParameters3D
{
public:
    /**
     * @brief Конструктор для изотропного тела
     * @param E Модуль Юнга
     * @param nu Коэффициент Пуассона
     */
    MechanicalParameters3D(double E, double nu);
    /**
     * @brief Конструктор для анизотропного тела (E, nu, G - независимые константы)
     * @param E Модуль Юнга
     * @param nu Коэффициент Пуассона
     * @param G Модуль сдвига
     */
    MechanicalParameters3D(double E, double nu, double G);
    /**
     * @brief Конструктор для анизотропного тела (ортотропия, общий случай, 9 независимых параметров)
     * @param E1 Модуль Юнга для первого направления (x)
     * @param E2 Модуль Юнга для второго направления (y)
     * @param E3 Модуль Юнга для третьего направления (z)
     * @param nu12 Коэффициент Пуассона в плоскости 1-2 (xy)
     * @param nu13 Коэффициент Пуассона в плоскости 1-3 (xz)
     * @param nu23 Коэффициент Пуассона в плоскости 2-3 (yz)
     * @param G12 Модуль сдвига в плоскости 1-2 (xy)
     * @param G13 Модуль сдвига в плоскости 1-3 (xz)
     * @param G23 Модуль сдвига в плоскости 2-3 (yz)
     */
    MechanicalParameters3D(double E1, double E2, double E3, double nu12, double nu13, double nu23, double G12, double G13, double G23);

    double E1() const;
    void setE1(double E1);

    double E2() const;
    void setE2(double E2);

    double E3() const;
    void setE3(double E3);

    double nu21() const;
    void setNu21(double nu21);

    double nu31() const;
    void setNu31(double nu31);

    double nu12() const;
    void setNu12(double nu12);

    double nu32() const;
    void setNu32(double nu32);

    double nu13() const;
    void setNu13(double nu13);

    double nu23() const;
    void setNu23(double nu23);

    double G23() const;
    void setG23(double G23);

    double G13() const;
    void setG13(double G13);

    double G12() const;
    void setG12(double G12);

private:
    double E1_;
    double E2_;
    double E3_;
    double nu21_;
    double nu31_;
    double nu12_;
    double nu32_;
    double nu13_;
    double nu23_;
    double G23_;
    double G13_;
    double G12_;
};

#endif // MECHANICALPARAMETERS3D_H
